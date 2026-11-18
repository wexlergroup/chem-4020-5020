/*
  Shared parts of the computational studios: the page frame, the colors of the book, and the styling of the charts.
  Every studio imports this file and kit.css, so that all studios have one look.

  Inside the notes, a studio follows the light or dark mode of the book and sets the height of its own frame.
  On paper it uses the colors of the light page.
  On its own page it follows the system, or ?theme=light or ?theme=dark.
  With ?static=1 it shows only the panels, for the still picture that the PDF of the notes uses.
*/
import React, { useEffect, useId, useRef, useState } from 'react';
import { flushSync } from 'react-dom';
import './kit.css';

const NAMES = ['ground', 'text', 'muted', 'border', 'fill', 'blue', 'vermilion', 'green', 'purple', 'gray', 'axis'];
const params = new URLSearchParams(window.location.search);
if (params.has('static')) document.body.classList.add('static');     // ?static=1: one still picture, for the PDF of the notes

function bookRoot() {                                 // the <html> element of the book, when the studio is in a frame
  try { return window.parent !== window ? window.parent.document.documentElement : null; } catch (e) { return null; }
}

let printing = false;                                 // on paper, the colors of the light page, whatever the mode of the screen

function readTheme() {
  const host = bookRoot();
  let dark = false;
  if (printing) dark = false;
  else if (params.get('theme')) dark = params.get('theme') === 'dark';
  else if (host) dark = host.classList.contains('dark');
  else { try { dark = window.matchMedia('(prefers-color-scheme: dark)').matches; } catch (e) { dark = false; } }
  document.documentElement.dataset.theme = dark ? 'dark' : 'light';
  const css = getComputedStyle(document.documentElement);
  const color = {};
  for (const name of NAMES) color[name] = css.getPropertyValue('--' + name).trim();
  return { dark, color };
}

function fit() {                                      // inside the notes, make the frame as tall as the content
  try {
    const frame = window.frameElement;
    if (!frame) return;
    const height = Math.ceil(document.documentElement.getBoundingClientRect().height) + 'px';
    frame.style.height = height;
    const box = frame.parentElement;                  // the book theme wraps a frame in a box with a fixed 5:3 shape
    if (box && box.style.paddingBottom) { box.style.paddingBottom = '0'; box.style.display = 'block'; }
    if (box && box.style.paddingBottom === '0px') box.style.height = height;
  } catch (e) { /* the studio is on its own page */ }
}

// A picture that is drawn in a loop of animation frames would be printed in its old colors. onPaper(redraw) calls
// redraw at once when the page is about to be printed, and again after it. It returns the function that ends this.
const REDRAW = 'studio-redraw';
export function onPaper(redraw) {
  window.addEventListener(REDRAW, redraw);
  return () => window.removeEventListener(REDRAW, redraw);
}

// On paper the frame keeps the width that it has on the screen, so that nothing in the studio is laid out again
// while the page is printed. The print rules of the notes (book/custom.css) read this width and reduce the frame
// to the width of the text.
function keepWidth(on) {
  try {
    const frame = window.frameElement;
    if (!frame) return;
    if (on) frame.style.setProperty('--paper-width', frame.getBoundingClientRect().width + 'px');
    else frame.style.removeProperty('--paper-width');
  } catch (e) { /* the studio is on its own page */ }
}

// The colors of the present mode, as { dark, color }. The component draws again when the mode changes.
export function useBook() {
  const [theme, setTheme] = useState(readTheme);
  useEffect(() => {
    try { if (window.frameElement) document.body.classList.add('embedded'); } catch (e) { /* own page */ }
    const update = () => setTheme(readTheme());
    const host = bookRoot();
    let observer = null, media = null;
    if (host) { observer = new MutationObserver(update); observer.observe(host, { attributes: true, attributeFilter: ['class'] }); }
    else { try { media = window.matchMedia('(prefers-color-scheme: dark)'); media.addEventListener('change', update); } catch (e) { media = null; } }
    let sizes = null;
    try { sizes = new ResizeObserver(fit); sizes.observe(document.body); } catch (e) { sizes = null; }
    window.addEventListener('load', fit);
    fit();
    const print = (on) => () => {                     // before and after the page is printed
      keepWidth(on);
      printing = on;
      const was = document.documentElement.dataset.theme === 'dark', next = readTheme();
      if (next.dark === was) return;                  // the light page is printed as it is
      flushSync(() => setTheme(next));
      window.dispatchEvent(new Event(REDRAW));
      document.body.getBoundingClientRect();          // lays the redrawn charts out now, or the printer leaves their curves out
    };
    const before = print(true), after = print(false);
    window.addEventListener('beforeprint', before);
    window.addEventListener('afterprint', after);
    return () => {
      if (observer) observer.disconnect();
      if (media) media.removeEventListener('change', update);
      if (sizes) sizes.disconnect();
      window.removeEventListener('load', fit);
      window.removeEventListener('beforeprint', before);
      window.removeEventListener('afterprint', after);
    };
  }, []);
  return theme;
}

// A variable, in italics: <V>T</V>, with a subscript: <V sub="B">k</V>, or a molar quantity with a bar: <V bar>v</V>
export function V({ children, sub, sup, bar }) {
  return <><i className={bar ? 'v bar' : 'v'}>{children}</i>{sub ? <sub>{sub}</sub> : null}{sup ? <sup>{sup}</sup> : null}</>;
}

export function Studio({ title, lede, children }) {
  return (
    <main>
      <h1>{title}</h1>
      <p className="lede">{lede}</p>
      {children}
    </main>
  );
}

export function Panel({ letter, title, children }) {
  return (
    <section className="panel">
      <h2><b>{letter}</b>{title}</h2>
      {children}
    </section>
  );
}

// What a screen reader says for a value on the page: the text as shown, with "10 to the power 3" for a power of ten.
function spoken(node) {
  let text = '';
  for (const part of node.childNodes) text += part.nodeName === 'SUP' && /10$/.test(text) ? ' to the power ' + part.textContent : part.textContent;
  return text;
}

// For a screen reader the slider is named by its label, and its value is read as it is shown, with the unit.
export function Slider({ label, value, show, min, max, step, onChange, disabled }) {
  const name = useId(), shown = useRef(null), input = useRef(null);
  useEffect(() => { input.current.setAttribute('aria-valuetext', spoken(shown.current)); });
  return (
    <label className={disabled ? 'off' : undefined}>
      <span id={name}>{label}</span> <output ref={shown}>{show}</output>
      <input ref={input} aria-labelledby={name} type="range" min={min} max={max} step={step} value={value} disabled={disabled} onChange={(e) => onChange(parseFloat(e.target.value))} />
    </label>
  );
}

export function Select({ label, value, options, onChange }) {          // options: [[value, text], ...]
  return (
    <label>
      {label}
      <select value={value} onChange={(e) => onChange(e.target.value)}>
        {options.map(([key, text]) => <option key={key} value={key}>{text}</option>)}
      </select>
    </label>
  );
}

export function Checks({ label, items, onChange }) {                    // items: [[key, text, checked], ...]
  return (
    <fieldset>
      <legend>{label}</legend>
      <div className="checks">
        {items.map(([key, text, checked]) => (
          <label key={key}><input type="checkbox" checked={checked} onChange={(e) => onChange(key, e.target.checked)} />{text}</label>
        ))}
      </div>
    </fieldset>
  );
}

// A number as m x 10^n with three significant figures: <Sci value={3.11e-19} />
export function Sci({ value, digits = 3 }) {
  if (value === 0) return <>0</>;
  let power = Math.floor(Math.log10(Math.abs(value)));
  let mantissa = Number((value / 10 ** power).toPrecision(digits));
  if (Math.abs(mantissa) >= 10) { mantissa /= 10; power += 1; }
  return <>{mantissa.toFixed(digits - 1)} × 10<sup>{String(power).replace('-', '−')}</sup></>;
}

// A random number generator with a seed, so that a studio starts in the same state every time
export function seeded(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export const STATIC = params.has('static');           // true for the still picture

// Styling of a recharts axis: open axes, no grid, ticks and labels in the muted color of the book
export function axis(color) {
  return {
    stroke: color.axis,
    tick: { fill: color.text, fontSize: 12 },
    tickLine: { stroke: color.axis },
    axisLine: { stroke: color.axis },
    tickMargin: 4,
  };
}
export function axisLabel(color, value, side) {
  if (side === 'left') return { value, angle: -90, position: 'insideLeft', offset: 12, fill: color.text, fontSize: 12.5, style: { textAnchor: 'middle' } };
  return { value, position: 'insideBottom', offset: -4, fill: color.text, fontSize: 12.5 };
}
export const aloud = (number) => String(number).replace(/^[-−]/, 'minus ');   // a number as a screen reader should say it
export const LINE = 2.4;                              // width of a data line
export const DASH = '4 3';                            // dashes of a reference line

// Tick values 1, 2, 5 times a power of ten between lo and hi, for a logarithmic axis
export function logTicks(lo, hi) {
  const ticks = [];
  for (let power = Math.floor(Math.log10(lo)); power <= Math.ceil(Math.log10(hi)); power++) {
    for (const m of [1, 2, 5]) {
      const v = Number((m * 10 ** power).toPrecision(1));
      if (v >= lo && v <= hi) ticks.push(v);
    }
  }
  return ticks;
}

// Round tick values from 0 to max, about n intervals, 1, 2, 2.5, or 5 times a power of ten apart
export function niceTicks(max, n) {
  const rough = max / n, power = 10 ** Math.floor(Math.log10(rough));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * power).find((s) => s >= rough - 1e-12);
  const ticks = [];
  for (let v = 0; v <= max + 1e-9; v += step) ticks.push(Math.round(v * 1e6) / 1e6);
  return ticks;
}
