(function(){
  function ensureToastContainer(){
    let c = document.getElementById('copy-toast-container');
    if (!c) {
      c = document.createElement('div');
      c.id = 'copy-toast-container';
      c.style.position = 'fixed';
      c.style.top = '12px';
      c.style.right = '12px';
      c.style.zIndex = '9999';
      c.style.pointerEvents = 'none';
      document.body.appendChild(c);
    }
    return c;
  }
  function showToast(msg){
    const c = ensureToastContainer();
    const t = document.createElement('div');
    t.textContent = msg;
    t.style.background = 'rgba(0,0,0,0.8)';
    t.style.color = '#fff';
    t.style.padding = '8px 12px';
    t.style.marginTop = '6px';
    t.style.borderRadius = '6px';
    t.style.fontFamily = 'system-ui, -apple-system, Segoe UI, Roboto, sans-serif';
    t.style.fontSize = '12px';
    t.style.boxShadow = '0 2px 8px rgba(0,0,0,0.3)';
    c.appendChild(t);
    setTimeout(()=>{ t.style.transition='opacity 250ms'; t.style.opacity='0'; }, 1200);
    setTimeout(()=>{ c.removeChild(t); }, 1600);
  }
  document.querySelectorAll('.js-plotly-plot').forEach(function(p){
    p.on('plotly_click', function(e){
      var geneId = e?.points?.[0]?.text;
      if (!geneId) return;
      const text = String(geneId);
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(function(){
          showToast('Copied Gene ID: ' + text);
        }).catch(function(){
          showToast('Copy failed');
        });
      } else {
        // Fallback for older browsers
        const ta = document.createElement('textarea');
        ta.value = text;
        document.body.appendChild(ta);
        ta.select();
        try { document.execCommand('copy'); showToast('Copied Gene ID: ' + text); }
        catch(err){ showToast('Copy failed'); }
        document.body.removeChild(ta);
      }
    });
  });
})();
