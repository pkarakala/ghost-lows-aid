/*
 * Test-only preload (node --require preload_clock.js): moves this process's wall clock.
 *
 *   GHOSTLOWS_AID_FAKE_NOW=<ISO time>   Date.now() and new Date() return that fixed instant
 *   GHOSTLOWS_AID_CLOCK_SHIFT_MS=<ms>   ... or the real time shifted by that many milliseconds
 *
 * Used to pin the shipped oref0 CLIs to the simulated clock (tests/test_bridge_cli.py) and to run
 * the bridge at shifted wall clocks (tests/test_bridge_clock.py). Never used by the simulation.
 */
'use strict';

const RealDate = Date;
const fixed = process.env.GHOSTLOWS_AID_FAKE_NOW
  ? RealDate.parse(process.env.GHOSTLOWS_AID_FAKE_NOW)
  : null;
const shift = Number(process.env.GHOSTLOWS_AID_CLOCK_SHIFT_MS || 0);
const now = () => (fixed !== null ? fixed : RealDate.now() + shift);

class ShiftedDate extends RealDate {
  constructor(...args) {
    if (args.length === 0) super(now());
    else super(...args);
  }
  static now() {
    return now();
  }
}
// Date() called without `new` returns the current time as a string.
global.Date = new Proxy(ShiftedDate, {
  apply() {
    return new RealDate(now()).toString();
  },
});
