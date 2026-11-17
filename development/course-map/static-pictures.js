// Draws the static pictures of the course map that stand in for it in the PDF of the notes:
// book/_static/course-map/overview.png (every module folded) and module-1.png to module-7.png.
// Run it again after any change to book/computational-studios/course-map.html:
//
//   npm install playwright && npx playwright install chromium   (once)
//   node development/course-map/static-pictures.js
//
const { chromium } = require("playwright");
const path = require("path");
const book = path.resolve(__dirname, "..", "..", "book");
const map = "file://" + path.join(book, "computational-studios", "course-map.html");
const out = path.join(book, "_static", "course-map");

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 700, height: 900 }, deviceScaleFactor: 2 });
  page.on("pageerror", error => { console.error(error.message); process.exitCode = 1; });
  for (let n = 0; n <= 7; n++) {                  // 0 is the map with every module folded
    await page.goto(`${map}?module=${n}&static=1&theme=light`);
    await page.locator(".roadmap").screenshot({ path: path.join(out, n ? `module-${n}.png` : "overview.png") });
  }
  await browser.close();
})();
