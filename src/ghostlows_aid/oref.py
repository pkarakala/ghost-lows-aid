"""The oref0 controller: a client for node/bridge.js and the controller step built on it.

The bridge is one long-lived node process per Python process, speaking JSON lines (TECHNICAL_DESIGN
§2, option B). It holds no state, so a crashed or hung bridge is restarted and the request resent.
Every request carries the simulated clock; the bridge never reads the wall clock (bridge.js header).
"""

from __future__ import annotations

import json
import os
import queue
import shutil
import subprocess
import threading
import time
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
NODE_DIR = Path(os.environ.get("GHOSTLOWS_AID_NODE_DIR") or REPO / "node")
BRIDGE_JS = NODE_DIR / "bridge.js"
# Never "fakecgm": determine-basal.js:177-185 skips the flat-CGM check for that device name.
DEVICE = "ghostlows-aid-cgm"
NO_TEMP = {"duration": 0, "rate": 0, "temp": "absolute"}  # what the pump reports with no temp


class BridgeError(RuntimeError):
    pass


def node_available() -> bool:
    return shutil.which("node") is not None and (NODE_DIR / "node_modules" / "oref0").is_dir()


class Bridge:
    """A node/bridge.js worker. `request` sends one JSON object and returns the reply."""

    def __init__(
        self,
        timeout_s: float = 60.0,
        preload: str | Path | None = None,
        env: dict[str, str] | None = None,
        retries: int = 1,
    ):
        self.timeout_s = timeout_s
        self.preload = preload  # a node --require script (tests shift the wall clock with it)
        self.env = env or {}
        self.retries = retries
        self.restarts = 0
        self._proc: subprocess.Popen | None = None
        self._lines: queue.Queue = queue.Queue()
        self._next_id = 0

    # -- process lifecycle --
    def start(self) -> None:
        node = shutil.which("node")
        if node is None:
            raise BridgeError("node is not on PATH")
        cmd = [node]
        if self.preload is not None:
            cmd += ["--require", str(self.preload)]
        cmd.append(str(BRIDGE_JS))
        env = {**os.environ, **self.env, "TZ": "UTC"}  # oref0 uses local-time getHours()
        self._proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            cwd=str(NODE_DIR),
            env=env,
            text=True,
            bufsize=1,
        )
        self._lines = queue.Queue()
        threading.Thread(target=self._read, args=(self._proc, self._lines), daemon=True).start()

    @staticmethod
    def _read(proc: subprocess.Popen, lines: queue.Queue) -> None:
        for line in proc.stdout:
            lines.put(line)
        lines.put(None)  # EOF: the process ended

    def close(self) -> None:
        proc, self._proc = self._proc, None
        if proc is None:
            return
        try:
            proc.stdin.write('{"op":"exit"}\n')
            proc.stdin.flush()
            proc.wait(timeout=5)
        except (OSError, subprocess.TimeoutExpired, ValueError):
            proc.kill()
            proc.wait()

    def __enter__(self) -> Bridge:
        self.start()
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __del__(self):  # pragma: no cover - best effort
        try:
            self.close()
        except Exception:
            pass

    # -- requests --
    def request(self, payload: dict) -> dict:
        attempts = 0
        while True:
            try:
                return self._request_once(payload)
            except BridgeError:
                attempts += 1
                if attempts > self.retries:
                    raise
                self._kill()
                self.restarts += 1

    def _kill(self) -> None:
        proc, self._proc = self._proc, None
        if proc is not None:
            proc.kill()
            proc.wait()

    def _request_once(self, payload: dict) -> dict:
        if self._proc is None or self._proc.poll() is not None:
            self.start()
        self._next_id += 1
        rid = self._next_id
        line = json.dumps({**payload, "id": rid}, separators=(",", ":"), allow_nan=False)
        try:
            self._proc.stdin.write(line + "\n")
            self._proc.stdin.flush()
        except (BrokenPipeError, OSError) as e:
            raise BridgeError(f"bridge stdin closed: {e}") from e
        deadline = time.monotonic() + self.timeout_s
        while True:
            try:
                raw = self._lines.get(timeout=max(0.0, deadline - time.monotonic()))
            except queue.Empty as e:
                raise BridgeError(f"bridge timed out after {self.timeout_s} s") from e
            if raw is None:
                raise BridgeError("bridge exited")
            reply = json.loads(raw)
            if reply.get("id") == rid:
                break
        if not reply.get("ok"):
            # an error inside oref0 is a result, not a crash: report it, do not retry
            raise OrefError(reply.get("error", "unknown bridge error"), reply.get("stack"))
        return reply


class OrefError(RuntimeError):
    """oref0 itself raised (bad inputs); the bridge is fine."""

    def __init__(self, message: str, stack: str | None = None):
        super().__init__(message)
        self.stack = stack


_shared: Bridge | None = None


def shared_bridge() -> Bridge:
    """The process-wide bridge, started on first use (one node worker per Python worker)."""
    global _shared
    if _shared is None:
        _shared = Bridge()
        _shared.start()
    return _shared


# ---- Nightscout-style records ----
def iso(ms: int) -> str:
    """Zoned ISO-8601 time with milliseconds, as Nightscout writes `created_at`/`dateString`."""
    return datetime.fromtimestamp(ms / 1000, tz=UTC).isoformat(timespec="milliseconds")[:-6] + "Z"


def glucose_entry(ms: int, value: float, device: str = DEVICE) -> dict:
    """One CGM reading as a Nightscout entry (lib/glucose-get-last.js:1-14 reads glucose or sgv,
    date and device; lib/determine-basal/cob.js:46-50 needs dateString)."""
    return {
        "date": int(ms),
        "dateString": iso(ms),
        "sgv": value,
        "glucose": value,
        "device": device,
        "type": "sgv",
    }


# ---- the controller ----
@dataclass
class Action:
    """What the pump is told: a new temp (rate U/h for duration min; duration 0 cancels), or
    nothing (rate None: let the current temp run), plus an optional SMB in units."""

    rate: float | None = None
    duration: int | None = None
    smb: float = 0.0


@dataclass
class OrefController:
    """oref0 0.7.1 through the bridge. Stateless apart from its settings: every call gets the
    whole observed history and ledger."""

    profile: dict
    microbolus: bool = False  # determine_basal's microBolusAllowed (oref1 SMB)
    # SOURCE: TECHNICAL_DESIGN §1.2 (autosens fixed to ratio 1 in v1)
    autosens: dict | None = field(default_factory=lambda: {"ratio": 1.0})
    glucose_hours: float = 24.0  # glucose history sent per call
    treatment_hours: float = 24.0  # pump history sent per call
    skip_meal: bool = False  # TECHNICAL_DESIGN §2.1 lever; off unless shown output-identical
    bridge: Bridge | None = None
    name: str = "oref0"

    def step(
        self, clock_ms: int, glucose: list[dict], treatments: list[dict], currenttemp: dict
    ) -> tuple[Action, dict]:
        """One determine-basal cycle. `glucose` and `treatments` are newest first, as oref0 wants
        them (lib/glucose-get-last.js:16, lib/iob/history.js:274-280)."""
        bridge = self.bridge or shared_bridge()
        reply = bridge.request(
            {
                "op": "step",
                "clock": int(clock_ms),
                "glucose": glucose,
                "treatments": treatments,
                "carbs": [],
                "profile": self.profile,
                "currenttemp": currenttemp,
                "autosens": self.autosens,
                "microBolusAllowed": self.microbolus,
                "skip_meal": self.skip_meal,
            }
        )
        rT = reply["rT"]
        if "error" in rT:
            raise OrefError(rT["error"])
        action = Action(
            rate=rT.get("rate"),
            duration=None if rT.get("duration") is None else int(round(rT["duration"])),
            smb=float(rT.get("units") or 0.0),
        )
        return action, reply
