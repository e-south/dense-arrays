/*
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/docs/assets/showcase/player.js

Play the published sequence movies within their fixed artwork coordinates.

Module Author(s): Eric J. South
--------------------------------------------------------------------------------
*/
(function () {
  function prepare(container) {
    const button = container.querySelector('button');
    const video = container.querySelector('video');
    if (!button || !video || button.dataset.ready) return;
    button.dataset.ready = 'true';
    button.addEventListener('click', async function () {
      if (button.dataset.fallback) {
        window.location.href = video.currentSrc || video.src;
        return;
      }
      if (!video.paused) {
        video.pause();
        button.textContent = 'Play';
        button.setAttribute('aria-pressed', 'false');
        return;
      }
      try {
        await video.play();
        video.hidden = false;
        button.textContent = 'Pause';
        button.setAttribute('aria-pressed', 'true');
      } catch (_) {
        button.textContent = 'Open video';
        button.dataset.fallback = 'true';
      }
    });
  }
  function initialize() {
    document.querySelectorAll('.sequence-playback').forEach(prepare);
    const selected = document.getElementById(window.location.hash.slice(1));
    if (selected && selected.matches('.sequence-playback')) selected.querySelector('button')?.click();
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', initialize);
  else initialize();
})();
