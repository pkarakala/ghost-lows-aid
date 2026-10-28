#!/usr/bin/env node
/*
 * Ghost Lows AID: the oref0 bridge.
 *
 * One long-lived, stateless worker. It reads one JSON request per line on stdin and writes one
 * JSON reply per line on stdout. It requires the library modules of the unmodified oref0@0.7.1
 * npm release directly from node_modules/oref0/lib; it never edits, wraps or monkey-patches them.
 * Python owns all state (glucose history, pump ledger, current temp); every request carries it.
 *
 * Run with TZ=UTC: oref0 uses local-time getHours()/getMinutes() (lib/profile/basal.js:21,
 * lib/profile/isf.js:13, lib/iob/history.js:14, determine-basal.js:925 via skip_neutral_temps).
 *
 * The shipped CLIs are bypassed on purpose: bin/oref0-determine-basal.js:105-114 silently turns
 * SMB off unless an `autotune` directory exists in the working directory. Each stage below is wired
 * exactly as its CLI wires it (file:line cited), including the CLI-only checks of oref0-meal.js.
 *
 * Wall clock. Every clock argument the libraries accept is set to the simulated clock:
 *   - lib/iob inputs.clock (lib/iob/index.js:39, lib/iob/history.js:175);
 *   - lib/meal inputs.clock (lib/meal/index.js:19);
 *   - determine_basal currentTime (determine-basal.js:133-148; without it deliverAt and the
 *     BG age use new Date()).
 * One read cannot be reached by an argument: lib/meal/total.js:22-31 builds the IOB inputs for
 * lib/determine-basal/cob.js without a clock, so cob.js:26 calls lib/iob/history.js with
 * inputs.clock undefined and history.js:175 `new Date(tz(undefined))` reads moment's clock, the
 * wall clock by default. Its value decides which records count as "in the future" (history.js:276)
 * and where the temp-cancel event goes (history.js:370-372), so COB depended on when the
 * simulation was run. The bridge sets moment.now (moment's documented clock hook; moment and
 * moment-timezone share one instance here) to the simulated clock for every request, which is
 * exactly what a live loop sees. All other `new Date()` reads in the reachable code are only taken
 * when a time argument is missing, and oref0 always passes one there (lib/iob/calculate.js:18 from
 * total.js:72; lib/profile/basal.js:11 from history.js:522 and cob.js:144; lib/profile/isf.js:10
 * from cob.js:126). tests/test_bridge_clock.py runs the same requests at two shifted wall clocks.
 *
 * Logging. oref0 writes dozens of lines per call with console.error and process.stderr.write.
 * They are dropped (or returned in `log` when a request sets debug: true); rT.reason is kept.
 */
'use strict';

const path = require('path');

const OREF = path.join(__dirname, 'node_modules', 'oref0');
const getLastGlucose = require(OREF + '/lib/glucose-get-last');
const generateIOB = require(OREF + '/lib/iob');
const generateMeal = require(OREF + '/lib/meal');
const determineBasal = require(OREF + '/lib/determine-basal/determine-basal');
const tempBasalFunctions = require(OREF + '/lib/basal-set-temp');
const generateProfile = require(OREF + '/lib/profile');
const sharedNodeUtils = require(OREF + '/bin/oref0-shared-node-utils');
// The moment instance oref0 itself resolves (lib/iob/history.js:3,6 and lib/meal/*.js:3).
const moment = require(require.resolve('moment', { paths: [OREF + '/lib/iob'] }));
const momentTz = require(require.resolve('moment-timezone', { paths: [OREF + '/lib/iob'] }));
if (moment !== momentTz) {
  throw new Error('moment-timezone does not share the moment instance oref0 uses');
}
const wallNow = moment.now;

// ---- logging: oref0's stderr chatter is dropped, or captured for debug requests ----
const stdoutWrite = process.stdout.write.bind(process.stdout);
let captured = null; // array while a debug request runs
function sink(...args) {
  if (captured !== null) {
    captured.push(args.map((a) => (typeof a === 'string' ? a : safeString(a))).join(' '));
  }
}
function safeString(a) {
  try {
    return JSON.stringify(a);
  } catch (e) {
    return String(a);
  }
}
console.error = sink;
console.log = sink; // nothing in oref0's lib writes to stdout; keep the protocol clean anyway
console.warn = sink;
console.info = sink;
process.stderr.write = (chunk) => {
  if (captured !== null) captured.push(String(chunk));
  return true;
};

// ---- helpers ----
function clone(x) {
  // Each stage gets fresh inputs, as each CLI reads its own JSON files: oref0 mutates inputs
  // (lib/iob/history.js:269-270 sets timestamp; lib/glucose-get-last.js:10,48; cob.js:143-144).
  return x === undefined ? undefined : JSON.parse(JSON.stringify(x));
}
function ms(t0) {
  return Number(process.hrtime.bigint() - t0) / 1e6;
}
function isoClock(clock) {
  // Zoned ISO string, as clock-zoned.json is (lib/iob/index.js:36 warns on an unzoned clock).
  return new Date(clock).toISOString();
}
function pinClock(req) {
  // moment.now is the only wall-clock read not reachable by an argument (see header).
  if (req.wallclock === true) {
    moment.now = wallNow; // negative control for tests only: leave moment on the wall clock
  } else {
    const clock = Number(req.clock);
    moment.now = () => clock;
  }
}

// ---- the stages, each wired as its CLI ----
function stageGlucoseStatus(glucose) {
  // bin/oref0-determine-basal.js:91,95: glucose_data = require(glucose.json); getLastGlucose(...)
  return getLastGlucose(clone(glucose));
}

function stageIOB(req) {
  // bin/oref0-calculate-iob.js:56-82: history, history24 (none), profile, clock, autosens
  const inputs = {
    history: clone(req.treatments),
    history24: null,
    profile: clone(req.profile),
    clock: isoClock(req.clock),
  };
  if (req.autosens) {
    inputs.autosens = clone(req.autosens); // calculate-iob.js:78-80
  }
  return generateIOB(inputs);
}

function stageMeal(req) {
  // bin/oref0-meal.js:30-137, pumphistory = the ledger, carbhistory = req.carbs ({} if absent,
  // oref0-meal.js:100). The CLI's own checks are kept: carb_ratio bounds (lines 88-92) and the
  // 36-reading minimum (lines 131-135).
  const profile = clone(req.profile);
  if (typeof profile.carb_ratio === 'undefined' || profile.carb_ratio < 3) {
    return { carbs: 0, mealCOB: 0, reason: 'carb_ratio ' + profile.carb_ratio + ' out of bounds' };
  }
  const glucose = clone(req.glucose);
  const inputs = {
    history: clone(req.treatments),
    profile: profile,
    basalprofile: clone(req.profile.basalprofile),
    clock: isoClock(req.clock),
    carbs: req.carbs === undefined ? {} : clone(req.carbs),
    glucose: glucose,
  };
  const meal = generateMeal(inputs);
  if (glucose.length < 36) {
    meal.mealCOB = 0;
    meal.reason = 'not enough glucose data to calculate carb absorption';
  }
  return meal;
}

function stageDetermineBasal(req, glucoseStatus, iobData, mealData) {
  // bin/oref0-determine-basal.js:91-94 read the inputs from JSON files, 142-145 autosens,
  // 165 reservoir (null), 180-184 currentTime, 213-215 the call itself.
  return determineBasal(
    glucoseStatus,
    clone(req.currenttemp),
    iobData,
    clone(req.profile),
    req.autosens ? clone(req.autosens) : null,
    mealData,
    tempBasalFunctions,
    req.microBolusAllowed === true,
    req.reservoir_data === undefined ? null : req.reservoir_data,
    new Date(Number(req.clock))
  );
}

// ---- ops ----
const ops = {
  exit() {
    return {};
  },

  ping(req) {
    return {
      versions: {
        oref0: require(OREF + '/package.json').version,
        moment: moment.version,
        node: process.version,
      },
      tz: process.env.TZ || null,
      utcOffsetMin: new Date(Number(req.clock || 0)).getTimezoneOffset(),
      wall_ms: Date.now(), // the bridge's own probe (tests move the wall clock); oref0 never sees it
    };
  },

  profile(req) {
    // bin/oref0-get-profile.js:217-244: preferences first, then the pump settings and schedules.
    const inputs = {};
    const preferences = req.preferences || {};
    for (const pref of Object.keys(preferences)) inputs[pref] = clone(preferences[pref]);
    inputs.max_iob = inputs.max_iob || 0;
    inputs.settings = clone(req.settings);
    inputs.targets = clone(req.targets);
    inputs.basals = clone(req.basals);
    inputs.isf = clone(req.isf);
    inputs.carbratio = clone(req.carbratio);
    inputs.temptargets = clone(req.temptargets || []);
    inputs.model = clone(req.model || {});
    inputs.autotune = {};
    const finalResult = sharedNodeUtils.initFinalResults();
    const profile = generateProfile(finalResult, inputs);
    if (profile === -1) throw new Error('lib/profile rejected the inputs: ' + finalResult.err);
    return { profile: profile };
  },

  glucose_status(req) {
    return { glucose_status: stageGlucoseStatus(req.glucose) };
  },

  iob(req) {
    return { iob: JSON.parse(JSON.stringify(stageIOB(req))) };
  },

  meal(req) {
    return { meal_data: JSON.parse(JSON.stringify(stageMeal(req))) };
  },

  determine_basal(req) {
    // Precomputed inputs, as the oref0 mocha fixtures pass them (used by the parity tests).
    const rT = determineBasal(
      clone(req.glucose_status),
      clone(req.currenttemp),
      clone(req.iob_data),
      clone(req.profile),
      req.autosens === undefined ? undefined : clone(req.autosens),
      clone(req.meal_data),
      tempBasalFunctions,
      req.microBolusAllowed,
      req.reservoir_data,
      req.currentTime === undefined ? undefined : new Date(req.currentTime)
    );
    return { rT: JSON.parse(JSON.stringify(rT)) };
  },

  step(req) {
    const timing = {};
    let t = process.hrtime.bigint();
    const glucoseStatus = stageGlucoseStatus(req.glucose);
    timing.glucose_status = ms(t);

    t = process.hrtime.bigint();
    // The CLI passes iob.json through a file: a JSON round trip keeps the same types (Dates as
    // ISO strings) that bin/oref0-determine-basal.js:93 reads back.
    const iobArray = JSON.parse(JSON.stringify(stageIOB(req)));
    timing.iob = ms(t);

    t = process.hrtime.bigint();
    let mealData;
    if (req.skip_meal === true) {
      // Lever under test (TECHNICAL_DESIGN §2.1); not used unless outputs are shown identical.
      mealData = { carbs: 0, nsCarbs: 0, bwCarbs: 0, journalCarbs: 0, mealCOB: 0,
        currentDeviation: 0, maxDeviation: 0, minDeviation: 0, slopeFromMaxDeviation: 0,
        slopeFromMinDeviation: 0, allDeviations: [], lastCarbTime: 0, bwFound: false };
    } else {
      mealData = JSON.parse(JSON.stringify(stageMeal(req)));
    }
    timing.meal = ms(t);

    t = process.hrtime.bigint();
    const rT = JSON.parse(JSON.stringify(stageDetermineBasal(req, glucoseStatus, iobArray, mealData)));
    timing.determine_basal = ms(t);

    const out = { glucose_status: glucoseStatus, iob: iobArray[0], meal_data: mealData, rT: rT };
    if (req.iob_array === true) out.iob_array = iobArray;
    out.timing_ms = timing;
    return out;
  },
};

function handle(line) {
  const t0 = process.hrtime.bigint();
  let req;
  try {
    req = JSON.parse(line);
  } catch (e) {
    return { ok: false, error: 'bad JSON request: ' + e.message };
  }
  const id = req.id === undefined ? null : req.id;
  const op = ops[req.op];
  if (!op) return { id: id, ok: false, error: 'unknown op ' + req.op };
  captured = req.debug === true ? [] : null;
  try {
    pinClock(req);
    const out = op(req);
    out.id = id;
    out.op = req.op;
    out.ok = true;
    if (captured !== null) out.log = captured;
    out.total_ms = ms(t0);
    return out;
  } catch (e) {
    return { id: id, ok: false, error: String(e && e.message ? e.message : e), stack: e && e.stack };
  } finally {
    captured = null;
    moment.now = wallNow;
  }
}

// ---- JSON lines over stdin/stdout ----
let buffer = '';
process.stdin.setEncoding('utf8');
process.stdin.on('data', (chunk) => {
  buffer += chunk;
  let nl;
  while ((nl = buffer.indexOf('\n')) >= 0) {
    const line = buffer.slice(0, nl);
    buffer = buffer.slice(nl + 1);
    if (line.trim() === '') continue;
    const reply = handle(line);
    stdoutWrite(JSON.stringify(reply) + '\n');
    if (reply.ok && reply.op === 'exit') process.exit(0);
  }
});
process.stdin.on('end', () => process.exit(0));
