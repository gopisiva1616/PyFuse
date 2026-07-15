(function () {
  function setSize(scale) {
    const clamped = Math.max(0.9, Math.min(1.4, scale));
    document.documentElement.style.setProperty('--pyfuse-font-scale', clamped.toString());
    try {
      localStorage.setItem('pyfuse-font-scale', clamped.toString());
    } catch (e) {
      // ignore storage errors in restricted browsers
    }
  }

  function getStoredScale() {
    try {
      const raw = localStorage.getItem('pyfuse-font-scale');
      if (!raw) {
        return 1.0;
      }
      const parsed = parseFloat(raw);
      if (Number.isNaN(parsed)) {
        return 1.0;
      }
      return parsed;
    } catch (e) {
      return 1.0;
    }
  }

  function buildControl() {
    if (document.querySelector('.pyfuse-font-controls')) {
      return;
    }

    const box = document.createElement('div');
    box.className = 'pyfuse-font-controls';

    const label = document.createElement('span');
    label.textContent = 'Font size';

    const smaller = document.createElement('button');
    smaller.type = 'button';
    smaller.textContent = 'A-';
    smaller.addEventListener('click', function () {
      setSize(getStoredScale() - 0.05);
    });

    const reset = document.createElement('button');
    reset.type = 'button';
    reset.textContent = 'A';
    reset.addEventListener('click', function () {
      setSize(1.0);
    });

    const larger = document.createElement('button');
    larger.type = 'button';
    larger.textContent = 'A+';
    larger.addEventListener('click', function () {
      setSize(getStoredScale() + 0.05);
    });

    box.appendChild(label);
    box.appendChild(smaller);
    box.appendChild(reset);
    box.appendChild(larger);
    document.body.appendChild(box);
  }

  document.addEventListener('DOMContentLoaded', function () {
    setSize(getStoredScale());
    buildControl();
  });
})();
