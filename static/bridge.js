let activeCharts = {};

document.addEventListener("DOMContentLoaded", () => {
    const fileInput = document.getElementById('fileInput');
    const processBtn = document.getElementById('processBtn');
    const downloadCsvBtn = document.getElementById('downloadCsvBtn');
    const downloadFaaBtn = document.getElementById('downloadFaaBtn');
    
    if (fileInput && processBtn) {
        fileInput.addEventListener('change', () => { 
            processBtn.disabled = !fileInput.files.length; 
        });
        processBtn.addEventListener('click', sendDataToPython);
    }

    // Hubungkan tombol unduh ke fungsi masing-masing
    if (downloadCsvBtn) downloadCsvBtn.addEventListener('click', exportCSV);
    if (downloadFaaBtn) downloadFaaBtn.addEventListener('click', exportFASTAProtein);

    // Otomatis memperbarui jika konfigurasi diubah setelah file terunggah
    ['lenFilter', 'winSize', 'frameSelector'].forEach(id => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('change', () => {
                if (fileInput && fileInput.files.length) sendDataToPython();
            });
        }
    });
    
    const lenFilter = document.getElementById('lenFilter');
    if (lenFilter) {
        lenFilter.addEventListener('input', (e) => {
            const el = document.getElementById('lenVal');
            if (el) el.innerText = e.target.value;
        });
    }
    
    const winSize = document.getElementById('winSize');
    if (winSize) {
        winSize.addEventListener('input', (e) => {
            const el = document.getElementById('winVal');
            if (el) el.innerText = e.target.value;
        });
    }
});

function sendDataToPython() {
    const fileInput = document.getElementById('fileInput');
    if (!fileInput || !fileInput.files.length) return;

    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('minLen', document.getElementById('lenFilter')?.value || 0);
    formData.append('winSize', document.getElementById('winSize')?.value || 100);
    formData.append('frame', document.getElementById('frameSelector')?.value || 0);

    fetch('/analyze', { method: 'POST', body: formData })
    .then(res => {
        if (!res.ok) throw new Error("Server Flask mengembalikan status " + res.status);
        return res.json();
    })
    .then(data => {
        if (data.error) return alert(data.error);
        updateUI(data);
    })
    .catch(err => {
        console.error("Gagal mengambil data:", err);
    });
}

function updateUI(data) {
    // Aktifkan tombol unduh secara aman
    const csvBtn = document.getElementById('downloadCsvBtn');
    const faaBtn = document.getElementById('downloadFaaBtn');
    if (csvBtn) csvBtn.disabled = false;
    if (faaBtn) faaBtn.disabled = false;

    // 1. Update Konten Teks Widget
    try {
        if (document.getElementById('statThermal')) document.getElementById('statThermal').innerText = data.summary.gc_percent + "% (Rasio GC)";
        if (document.getElementById('statCodons')) document.getElementById('statCodons').innerText = `${data.summary.start_found} / ${data.summary.stop_found}`;
        if (document.getElementById('statN50')) document.getElementById('statN50').innerText = data.summary.total_bp + " bp";
    } catch(e) { console.error(e); }

    // 2. Update Viewer Sekuens
    try {
        if (document.getElementById('dnaView')) document.getElementById('dnaView').innerText = data.viewer.dna;
        if (document.getElementById('mrnaView')) document.getElementById('mrnaView').innerText = data.viewer.mrna;
        if (document.getElementById('peptideView')) document.getElementById('peptideView').innerText = data.viewer.peptide;
    } catch(e) { console.error(e); }

    // 3. Render Tabel Hasil Analisis
    try {
        const tbody = document.querySelector('#mainTable tbody');
        if (tbody) {
            tbody.innerHTML = '';
            data.table.forEach(row => {
                tbody.innerHTML += `<tr>
                    <td><strong>${row.id}</strong></td>
                    <td>${row.length}</td>
                    <td>${row.gc}%</td>
                    <td>${row.at}%</td>
                    <td>${row.skew}</td>
                    <td>${row.entropy}</td>
                    <td>${row.orf_len}</td>
                </tr>`;
            });
        }
    } catch(e) { console.error(e); }

    // 4. Render Grafik (Diletakkan paling akhir agar tidak mengunci fungsi lain jika gagal)
    try {
        renderCharts(data.charts);
    } catch(e) { console.error("Gagal memproses Chart.js:", e); }
}

function renderCharts(chartsData) {
    // Menghancurkan instansi chart lama secara aman sebelum membuat yang baru
    ['skewChart', 'purineChart', 'hydroChart', 'aminoChart'].forEach(id => {
        if (activeCharts[id]) {
            try { activeCharts[id].destroy(); } catch(e) {}
        }
    });

    if (typeof Chart === 'undefined') {
        console.error("Pustaka Chart.js tidak terdeteksi! Pastikan komputer terkoneksi internet.");
        return;
    }

    try {
        const ctxSkew = document.getElementById('skewChart');
        if (ctxSkew) {
            activeCharts['skewChart'] = new Chart(ctxSkew.getContext('2d'), {
                type: 'line',
                data: { labels: chartsData.skew_labels, datasets: [{ label: 'GC Skew', data: chartsData.skew_data, borderColor: '#6366f1', tension: 0.2 }] },
                options: { responsive: true, maintainAspectRatio: false }
            });
        }
    } catch(e) { console.error(e); }

    try {
        const ctxPurine = document.getElementById('purineChart');
        if (ctxPurine) {
            activeCharts['purineChart'] = new Chart(ctxPurine.getContext('2d'), {
                type: 'bar',
                data: { labels: ['Purin (A+G)', 'Pirimidin (C+T)'], datasets: [{ label: 'Basa Nitrogen', data: chartsData.purines, backgroundColor: ['#10b981', '#f43f5e'] }] },
                options: { responsive: true, maintainAspectRatio: false }
            });
        }
    } catch(e) { console.error(e); }

    try {
        const ctxHydro = document.getElementById('hydroChart');
        if (ctxHydro) {
            activeCharts['hydroChart'] = new Chart(ctxHydro.getContext('2d'), {
                type: 'line',
                data: { labels: Array.from({length: chartsData.hydro_data.length}, (_, i) => i + 1), datasets: [{ label: 'Skor Hidrofobisitas', data: chartsData.hydro_data, backgroundColor: 'rgba(245, 158, 11, 0.2)', borderColor: '#f59e0b', fill: true }] },
                options: { responsive: true, maintainAspectRatio: false }
            });
        }
    } catch(e) { console.error(e); }

    try {
        const ctxAmino = document.getElementById('aminoChart');
        if (ctxAmino) {
            activeCharts['aminoChart'] = new Chart(ctxAmino.getContext('2d'), {
                type: 'radar',
                data: { labels: ['Polar', 'Non-Polar', 'Kodon AUG', 'Kodon Stop'], datasets: [{ label: 'Proporsi Sifat Asam Amino', data: chartsData.amino_prop, borderColor: '#38bdf8', backgroundColor: 'rgba(56, 189, 248, 0.2)' }] },
                options: { responsive: true, maintainAspectRatio: false }
            });
        }
    } catch(e) { console.error(e); }
}

// === SISTEM EKSPOR / UNDUH DATA ===
function exportCSV() {
    const table = document.getElementById('mainTable');
    if (!table) return;
    
    let csv = "ID Sekuens,Panjang (bp),GC (%),AT (%),GC Skew,Entropy,ORF Length (aa)\n";
    const rows = table.querySelectorAll('tbody tr');
    
    rows.forEach(row => {
        const cols = row.querySelectorAll('td');
        if (cols.length > 1) {
            let rowData = Array.from(cols).map(col => col.innerText.replace('%', '')).join(',');
            csv += rowData + "\n";
        }
    });
    
    triggerDownload("data:text/csv;charset=utf-8," + encodeURIComponent(csv), 'laporan_sekuens.csv');
}

function exportFASTAProtein() {
    const peptideText = document.getElementById('peptideView')?.innerText;
    if (!peptideText || peptideText === "...") {
        alert("Belum ada data sekuens protein untuk diunduh.");
        return;
    }
    
    let fasta = `>Hasil_Translasi_Protein\n${peptideText}\n`;
    triggerDownload("data:text/plain;charset=utf-8," + encodeURIComponent(fasta), 'sekuens_protein.faa');
}

function triggerDownload(uri, filename) {
    let link = document.createElement('a'); 
    link.setAttribute('href', uri);
    link.setAttribute('download', filename); 
    document.body.appendChild(link);
    link.click(); 
    document.body.removeChild(link);
}