#!/usr/bin/env node
/*
 * Parity fixtures for the oref0 bridge (tests/test_bridge_parity.py).
 *
 * 1. Runs oref0's own tests/determine-basal.test.js in this process with stand-in describe/it,
 *    and records the arguments of every determine_basal call it makes. The recording wrapper is
 *    placed in this process's require cache only; no oref0 file is read differently or changed.
 * 2. Adds derived fixtures for the branches that file does not reach (SMB, stale data, flat CGM,
 *    the "fakecgm" simulator bypass, high noise), built from its standard initial conditions and a
 *    real lib/iob array.
 * 3. Calls the real determine_basal directly on each fixture with currentTime pinned (otherwise
 *    deliverAt is the wall clock, determine-basal.js:133-136) and prints
 *    {currentTime, fixtures: [{name, args, rT}]} as JSON on stdout.
 *
 * Usage: TZ=UTC node parity_fixtures.js > fixtures.json
 */
'use strict';

const path = require('path');

const OREF = path.join(__dirname, 'node_modules', 'oref0');
const DB = require.resolve(OREF + '/lib/determine-basal/determine-basal');
const realDetermineBasal = require(DB);
const tempBasalFunctions = require(OREF + '/lib/basal-set-temp');
const generateIOB = require(OREF + '/lib/iob');
const moment = require(require.resolve('moment', { paths: [OREF + '/lib/iob'] }));

const CURRENT_TIME = '2026-01-05T12:00:00.000Z';
const NOW = Date.parse(CURRENT_TIME);
const out = console.log.bind(console);
console.error = () => {};
console.log = () => {};
process.stderr.write = () => true;

const NAMES = ['glucose_status', 'currenttemp', 'iob_data', 'profile', 'autosens', 'meal_data',
  'tempBasalFunctions', 'microBolusAllowed', 'reservoir_data', 'currentTime'];

function snapshot(args) {
  const rec = {};
  args.forEach((a, i) => {
    if (NAMES[i] === 'tempBasalFunctions' || a === undefined) return;
    rec[NAMES[i]] = JSON.parse(JSON.stringify(a));
  });
  return rec;
}

// ---- 1. record the calls of oref0's own fixture file ----
const recorded = [];
let currentTest = null;
require.cache[DB].exports = function recordingDetermineBasal(...args) {
  recorded.push({ name: currentTest, args: snapshot(args) });
  return realDetermineBasal(...args);
};
global.describe = (name, fn) => fn();
global.it = (name, fn) => {
  currentTest = name;
  try {
    fn();
  } catch (e) {
    // oref0's own assertions: their pass/fail is mocha's business (84 passing), not ours
  }
  currentTest = null;
};
require(OREF + '/tests/determine-basal.test.js');
require.cache[DB].exports = realDetermineBasal;

// ---- 2. derived fixtures for the SMB, stale, flat, fakecgm and noise branches ----
// Standard initial conditions of tests/determine-basal.test.js:60-66, extended with the SMB
// preferences of lib/profile/index.js defaults.
const baseProfile = {
  max_iob: 2.5, dia: 3, type: 'current', current_basal: 0.9, max_daily_basal: 1.3, max_basal: 3.5,
  max_bg: 120, min_bg: 110, sens: 40, carb_ratio: 10,
  maxSMBBasalMinutes: 30, maxUAMSMBBasalMinutes: 30, SMBInterval: 3, bolus_increment: 0.1,
  enableUAM: true, enableSMB_with_COB: true, curve: 'rapid-acting', maxCOB: 120,
  min_5m_carbimpact: 8, remainingCarbsCap: 90, remainingCarbsFraction: 1,
  basalprofile: [{ i: 0, start: '00:00:00', minutes: 0, rate: 0.9 }],
};
const iso = (minAgo) => new Date(NOW - minAgo * 60000).toISOString();
const history = [
  { eventType: 'Temp Basal', rate: 1.5, duration: 30, created_at: iso(10), enteredBy: 'parity' },
  { eventType: 'Meal Bolus', insulin: 2.0, carbs: 40, created_at: iso(40), enteredBy: 'parity' },
];
moment.now = () => NOW;
const iobArray = JSON.parse(JSON.stringify(generateIOB({
  history: JSON.parse(JSON.stringify(history)), history24: null,
  profile: JSON.parse(JSON.stringify(baseProfile)), clock: CURRENT_TIME,
})));
const runningHigh = { duration: 20, rate: 1.5, temp: 'absolute' };
const meal = { carbs: 40, nsCarbs: 40, bwCarbs: 0, journalCarbs: 0, mealCOB: 30,
  currentDeviation: 4, maxDeviation: 6, minDeviation: 1, slopeFromMaxDeviation: -0.5,
  slopeFromMinDeviation: 0.8, allDeviations: [4, 5, 6, 3], lastCarbTime: NOW - 40 * 60000,
  bwFound: false };
const gs = (glucose, delta, shortAvg, longAvg, minAgo, extra) => Object.assign({
  delta: delta, glucose: glucose, short_avgdelta: shortAvg, long_avgdelta: longAvg,
  date: NOW - minAgo * 60000, noise: null, last_cal: 0, device: 'share2' }, extra || {});
const derived = [
  ['SMB: rising with COB, microbolus', gs(160, 6, 5, 4, 1), runningHigh, true, meal],
  ['SMB: no microbolus when microBolusAllowed is false', gs(160, 6, 5, 4, 1), runningHigh, false, meal],
  ['SMB: steep rise above maxDelta_bg_threshold disables SMB', gs(120, 30, 20, 10, 1), runningHigh, true, meal],
  ['stale: 20 min old BG with a high temp running', gs(160, 3, 3, 3, 20), runningHigh, true, meal],
  ['stale: 20 min old BG with a long zero temp', gs(80, -3, -3, -3, 20), { duration: 90, rate: 0, temp: 'absolute' }, false, meal],
  ['stale: 20 min old BG, no temp running', gs(160, 3, 3, 3, 20), { duration: 0, rate: 0, temp: 'absolute' }, false, meal],
  ['flat: unchanged CGM with a high temp running', gs(150, 0, 0.4, -0.6, 1), runningHigh, false, meal],
  ['flat: unchanged CGM, device fakecgm bypasses the check', gs(150, 0, 0.4, -0.6, 1, { device: 'fakecgm' }), runningHigh, false, meal],
  ['noise: noise 3 with a high temp running', gs(160, 3, 3, 3, 1, { noise: 3 }), runningHigh, false, meal],
  ['low suspend: low and falling with an IOB array', gs(72, -4, -4, -3, 1), runningHigh, true, meal],
];
for (const [name, glucoseStatus, currenttemp, smb, mealData] of derived) {
  recorded.push({ name: 'derived: ' + name, args: {
    glucose_status: glucoseStatus, currenttemp: currenttemp, iob_data: iobArray,
    profile: baseProfile, autosens: { ratio: 1.0 }, meal_data: mealData, microBolusAllowed: smb,
    reservoir_data: null } });
}

// ---- 3. the direct call, with currentTime pinned ----
const clone = (x) => (x === undefined ? undefined : JSON.parse(JSON.stringify(x)));
const fixtures = recorded.map((f) => {
  const a = f.args;
  const rT = realDetermineBasal(clone(a.glucose_status), clone(a.currenttemp), clone(a.iob_data),
    clone(a.profile), clone(a.autosens), clone(a.meal_data), tempBasalFunctions,
    a.microBolusAllowed, a.reservoir_data === undefined ? undefined : a.reservoir_data,
    new Date(CURRENT_TIME));
  return { name: f.name, args: a, rT: JSON.parse(JSON.stringify(rT)) };
});
out(JSON.stringify({ currentTime: CURRENT_TIME, fixtures: fixtures }));
