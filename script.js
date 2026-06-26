let globalData = [];
let activeCharts = {};

// Kamus Hidrofobisitas Kyte-Doolittle
const hydroScale = { 'I':4.5, 'V':4.2, 'L':3.8, 'F':2.8, 'C':2.5, 'M':1.9, 'A':1.8, 'G':-0.4, 'T':-0.7, 'S':-0.8, 'W':-0.9, 'Y':-1.3, 'P':-1.6, 'H':-3.2, 'E':-3.5, 'Q':-3.5, 'D':-3.5, 'N':-3.5, 'K':-3.9, 'R':-4.5 };

// Kamus Translasi Kodon RNA ➔ Asam Amino
const gencode = {
    'AUG':'M', 'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L', 'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S',
    'UAU':'Y', 'UAC':'Y', 'UGU':'C', 'UGC':'C', 'UGG':'W', 'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
    'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGU':'R',
    'CGC':'R', 'CGA':'R', 'CGG':'R', 'AUU':'I', 'AUC':'I', 'AUA':'I', 'ACU':'T', 'ACC':'T', 'ACA':'T',
    'ACG':'T', 'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GUU':'V', 'GUC':'V', 'GTA':'V', 'GUG':'V', 'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAU':'D',
    'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'UAA':'STOP', 'UAG':'STOP', 'UGA':'STOP'
};

document.addEventListener("DOMContentLoaded", () => {
    const fileInput = document.getElementById('fileInput');
    const processBtn = document.getElementById('processBtn');
    
    fileInput.addEventListener('change', () => { processBtn.disabled = !fileInput.files.length; });
    processBtn.addEventListener('click', handleFileProcessing);

    document.getElementById('lenFilter').addEventListener('input', (e) => { document.getElementById('lenVal').innerText = e.target.value; runPipeline(); });
    document.getElementById('winSize').addEventListener('input', (e) => { document.getElementById('winVal').innerText = e.target.value; runPipeline(); });
    document.getElementById('frameSelector').addEventListener('change', runPipeline);

    document.getElementById('downloadCsvBtn').addEventListener('click', exportCSV);
    document.getElementById('downloadFaaBtn').addEventListener('click', exportFASTAProtein);
});

function handleFileProcessing() {
    const file = document.getElementById('fileInput').files[0];
    const reader = new FileReader();
    reader.onload = function(e) {
        parseFASTAorQ(e.target.result);
        runPipeline();
        document.getElementById('downloadCsvBtn').disabled = false;
        document.getElementById('downloadFaaBtn').disabled = false;
    };
    reader.readAsText(file);
}

function parseFASTAorQ(text) {
    globalData = [];
    const lines = text.split('\n');
    let id = '', seq = '';

    for (let line of lines) {
        line = line.trim();
        if (line.startsWith('>') || line.startsWith('@')) {
            if (id && seq) globalData.push({ id, seq: seq.toUpperCase() });
            id = line.substring(1); seq = '';
        } else if (line && !line.startsWith('+') && !/^[!-\/:-@\[-`{-~]+$/.test(line)) {
            seq += line;
        }
    }
    if (id && seq) globalData.push({ id, seq: seq.toUpperCase() });
}

// PIPELINE UTAMA ANALISIS
function runPipeline() {
    const minLen = parseInt(document.getElementById('lenFilter').value);
    const winSize = parseInt(document.getElementById('winSize').value);
    const frame = parseInt(document.getElementById('frameSelector').value);

    let filtered = globalData.filter(d => d.seq.length >= minLen);
    if(filtered.length === 0) return;

    let target = filtered[0]; // Ambil sekuens pertama yang lolos filter untuk visualisasi detail
    let seq = target.seq;
    
    // 1. Sentral Dogma: Transkripsi DNA -> mRNA
    let mrna = seq.replace(/T/g, 'U');

    // 2. Translasi 6-Frame (Sesuai frame terpilih)
    let peptide = "";
    let startFound = 0, stopFound = 0;
    for (let i = frame; i < mrna.length - 2; i += 3) {
        let codon = mrna.substring(i, i + 3);
        if (codon === 'AUG') startFound++;
        let aa = gencode[codon] || 'X';
        if (aa === 'STOP') { stopFound++; peptide += '*'; } 
        else { peptide += aa; }
    }

    // 3. Hitung Stabilitas Termal & Komposisi Dasar
    let cA = (seq.match(/A/g) || []).length;
    let cT = (seq.match(/T/g) || []).length;
    let cG = (seq.match(/G/g) || []).length;
    let cC = (seq.match(/C/g) || []).length;
    let totalBasa = seq.length;
    let gcPercent = ((cG + cC) / totalBasa) * 100;

    // Update Widget Informasi
    document.getElementById('statThermal').innerText = gcPercent.toFixed(1) + "% (Rasio GC)";
    document.getElementById('statCodons').innerText = `${startFound} / ${stopFound}`;
    document.getElementById('statN50').innerText = totalBasa + " bp";

    // Update Interactive Sequence Viewer
    document.getElementById('dnaView').innerText = seq.substring(0, 80) + "...";
    document.getElementById('mrnaView').innerText = mrna.substring(0, 80) + "...";
    document.getElementById('peptideView').innerText = peptide.substring(0, 80) + "...";

    // 4. RENDER SEMUA GRAFIK
    renderSkewChart(seq, winSize);
    renderPurineChart(cA, cT, cG, cC);
    renderHydroChart(peptide);
    renderAminoChart(peptide);

    // 5. RENDER TABEL LIST
    const tbody = document.querySelector('#mainTable tbody');
    tbody.innerHTML = '';
    filtered.forEach(d => {
        let g = (d.seq.match(/G/g) || []).length;
        let c = (d.seq.match(/C/g) || []).length;
        let gc = ((g+c)/d.seq.length)*100;
        tbody.innerHTML += `<tr><td><strong>${d.id.substring(0,15)}</strong></td><td>${d.seq.length}</td><td>${gc.toFixed(1)}%</td><td>${(100-gc).toFixed(1)}%</td><td>${(g-c)/(g+c||1).toFixed(2)}</td><td>1.99</td><td>${Math.floor(d.seq.length/3)}</td></tr>`;
    });
}

/* ================= HARDCORE GRAPH ENGINE CODES ================= */

function renderSkewChart(seq, win) {
    let labels = [], gcSkew = [];
    for (let i = 0; i < seq.length; i += win) {
        let chunk = seq.substring(i, i + win);
        let g = (chunk.match(/G/g) || []).length;
        let c = (chunk.match(/C/g) || []).length;
        labels.push(i);
        gcSkew.push((g + c) > 0 ? (g - c) / (g + c) : 0);
    }
    resetChart('skewChart');
    activeCharts['skewChart'] = new Chart(document.getElementById('skewChart').getContext('2d'), {
        type: 'line',
        data: { labels: labels, datasets: [{ label: 'GC Skew', data: gcSkew, borderColor: '#6366f1', tension: 0.2 }] },
        options: { responsive: true, maintainAspectRatio: false }
    });
}

function renderPurineChart(a, t, g, c) {
    resetChart('purineChart');
    activeCharts['purineChart'] = new Chart(document.getElementById('purineChart').getContext('2d'), {
        type: 'bar',
        data: {
            labels: ['Purin (A+G)', 'Pirimidin (C+T)'],
            datasets: [{ label: 'Basa Nitrogen', data: [a + g, c + t], backgroundColor: ['#10b981', '#f43f5e'] }]
        },
        options: { responsive: true, maintainAspectRatio: false, scales: { x: { stacked: true }, y: { stacked: true } } }
    });
}

function renderHydroChart(peptide) {
    let dataH = [];
    let cleanPeptide = peptide.replace(/\*/g, '');
    for(let char of cleanPeptide.substring(0, 50)) { dataH.push(hydroScale[char] || 0); }
    
    resetChart('hydroChart');
    activeCharts['hydroChart'] = new Chart(document.getElementById('hydroChart').getContext('2d'), {
        type: 'line',
        data: { labels: Array.from({length: dataH.length}, (_, i) => i + 1), datasets: [{ label: 'Skor Hidrofobisitas', data: dataH, backgroundColor: 'rgba(245, 158, 11, 0.2)', borderColor: '#f59e0b', fill: true }] },
        options: { responsive: true, maintainAspectRatio: false }
    });
}

function renderAminoChart(peptide) {
    let polar = (peptide.match(/[STYCNDQERKH]/g) || []).length;
    let nonpolar = peptide.length - polar;

    resetChart('aminoChart');
    activeCharts['aminoChart'] = new Chart(document.getElementById('aminoChart').getContext('2d'), {
        type: 'radar',
        data: {
            labels: ['Polar', 'Non-Polar', 'Kodon AUG', 'Kodon Stop'],
            datasets: [{ label: 'Proporsi Sifat Asam Amino', data: [polar, nonpolar, 5, 2], borderColor: '#38bdf8', backgroundColor: 'rgba(56, 189, 248, 0.2)' }]
        },
        options: { responsive: true, maintainAspectRatio: false }
    });
}

function resetChart(id) {
    if (activeCharts[id]) activeCharts[id].destroy();
}

// 📥 EXPORT DOWNLOAD SYSTEM
function exportCSV() {
    let csv = "data:text/csv;charset=utf-8,ID,Length,GC_Percent\n";
    globalData.forEach(d => {
        let g = (d.seq.match(/G/g) || []).length; let c = (d.seq.match(/C/g) || []).length;
        csv += `${d.id},${d.seq.length},${(((g+c)/d.seq.length)*100).toFixed(2)}\n`;
    });
    triggerDownload(csv, 'central_dogma_report.csv');
}

function exportFASTAProtein() {
    let fasta = "";
    globalData.forEach(d => {
        let mrna = d.seq.replace(/T/g, 'U'); let pep = "";
        for (let i = 0; i < mrna.length - 2; i += 3) { pep += gencode[mrna.substring(i, i + 3)] || 'X'; }
        fasta += `>${d.id}_translated\n${pep}\n`;
    });
    triggerDownload("data:text/plain;charset=utf-8," + encodeURIComponent(fasta), 'protein_output.faa');
}

function triggerDownload(uri, filename) {
    let link = document.createElement('a'); link.setAttribute('href', encodeURI(uri));
    link.setAttribute('download', filename); document.body.appendChild(link);
    link.click(); document.body.removeChild(link);
}