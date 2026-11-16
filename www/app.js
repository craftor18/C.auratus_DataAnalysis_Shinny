(function () {
  var reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function countUp(el) {
    if (el.dataset.done) return;
    el.dataset.done = "1";
    var target = parseFloat(el.dataset.target) || 0;
    var fmt = function (v) { return Math.round(v).toLocaleString("en-US"); };
    if (reduceMotion) { el.textContent = fmt(target); return; }
    var start = null, dur = 1400;
    function step(ts) {
      if (!start) start = ts;
      var t = Math.min((ts - start) / dur, 1);
      var eased = 1 - Math.pow(1 - t, 3);
      el.textContent = fmt(target * eased);
      if (t < 1) requestAnimationFrame(step);
    }
    requestAnimationFrame(step);
  }

  var observer = "IntersectionObserver" in window ? new IntersectionObserver(function (entries) {
    entries.forEach(function (entry) {
      if (!entry.isIntersecting) return;
      entry.target.classList.add("is-visible");
      entry.target.querySelectorAll(".count-up").forEach(countUp);
      observer.unobserve(entry.target);
    });
  }, { threshold: 0.15 }) : null;

  function watch(root) {
    (root || document).querySelectorAll(".reveal:not(.is-visible)").forEach(function (el) {
      if (observer) observer.observe(el);
      else { el.classList.add("is-visible"); el.querySelectorAll(".count-up").forEach(countUp); }
    });
  }

  document.addEventListener("DOMContentLoaded", function () { watch(document); });

  document.addEventListener("click", function (e) {
    var link = e.target.closest("[data-nav]");
    if (!link || !window.Shiny) return;
    e.preventDefault();
    Shiny.setInputValue("goto_tab", link.dataset.nav, { priority: "event" });
    window.scrollTo({ top: 0, behavior: reduceMotion ? "auto" : "smooth" });
  });

  document.addEventListener("shown.bs.tab", function () { watch(document); });
})();
