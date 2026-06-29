from flask import Flask, render_template, request, jsonify
import re
import math

app = Flask(__name__)

# Kamus Hidrofobisitas Kyte-Doolittle
HYDRO_SCALE = { 'I':4.5, 'V':4.2, 'L':3.8, 'F':2.8, 'C':2.5, 'M':1.9, 'A':1.8, 'G':-0.4, 'T':-0.7, 'S':-0.8, 'W':-0.9, 'Y':-1.3, 'P':-1.6, 'H':-3.2, 'E':-3.5, 'Q':-3.5, 'D':-3.5, 'N':-3.5, 'K':-3.9, 'R':-4.5 }

# Kamus Translasi Kodon RNA
GENCODE = {
    'AUG':'M', 'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L', 'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S',
    'UAU':'Y', 'UAC':'Y', 'UGU':'C', 'UGC':'C', 'UGG':'W', 'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
    'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGU':'R',
    'CGC':'R', 'CGA':'R', 'CGG':'R', 'AUU':'I', 'AUC':'I', 'AUA':'I', 'ACU':'T', 'ACC':'T', 'ACA':'T',
    'ACG':'T', 'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GUU':'V', 'GUC':'V', 'GTA':'V', 'GUG':'V', 'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAU':'D',
    'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'UAA':'STOP', 'UAG':'STOP', 'UGA':'STOP'
}

def parse_fasta_or_q(text):
    sequences = []
    lines = text.split('\n')
    current_id = ''
    current_seq = []

    for line in lines:
        line = line.strip()
        if line.startswith('>') or line.startswith('@'):
            if current_id and current_seq:
                sequences.append({'id': current_id, 'seq': "".join(current_seq).upper()})
            current_id = line[1:]
            current_seq = []
        elif line and not line.startswith('+') and not re.match(r'^[!-\/:-@\[-`{-~]+$', line):
            current_seq.append(line)
            
    if current_id and current_seq:
        sequences.append({'id': current_id, 'seq': "".join(current_seq).upper()})
    return sequences

def hitung_entropi(seq):
    total = len(seq)
    if total == 0: return 0
    counts = {basa: seq.count(basa) for basa in set(seq)}
    entropy = 0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return round(entropy, 2)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze():
    file = request.files.get('file')
    min_len = int(request.form.get('minLen', 0))
    win_size = int(request.form.get('winSize', 100))
    frame = int(request.form.get('frame', 0))

    if not file:
        return jsonify({'error': 'No file uploaded'}), 400

    text = file.read().decode('utf-8')
    all_data = parse_fasta_or_q(text)
    
    filtered_data = [d for d in all_data if len(d['seq']) >= min_len]
    if not filtered_data:
        return jsonify({'error': 'No sequences match the criteria'}), 400

    target = filtered_data[0]
    seq = target['seq']
    
    mrna = seq.replace('T', 'U')

    peptide = ""
    start_found = 0
    stop_found = 0
    for i in range(frame, len(mrna) - 2, 3):
        codon = mrna[i:i+3]
        if codon == 'AUG': 
            start_found += 1
        aa = GENCODE.get(codon, 'X')
        if aa == 'STOP':
            stop_found += 1
            peptide += '*'
        else:
            peptide += aa

    g_count = seq.count('G')
    c_count = seq.count('C')
    a_count = seq.count('A')
    t_count = seq.count('T')
    total_basa = len(seq)
    gc_percent = ((g_count + c_count) / total_basa) * 100 if total_basa > 0 else 0

    skew_labels = []
    gc_skew_data = []
    for i in range(0, len(seq), win_size):
        chunk = seq[i:i+win_size]
        g = chunk.count('G')
        c = chunk.count('C')
        skew_labels.append(i)
        gc_skew_data.append((g - c) / (g + c) if (g + c) > 0 else 0)

    clean_peptide = peptide.replace('*', '')[:50]
    hydro_data = [HYDRO_SCALE.get(aa, 0) for aa in clean_peptide]

    polar_count = len(re.findall(r'[STYCNDQERKH]', peptide))
    nonpolar_count = len(peptide) - polar_count

    table_rows = []
    for d in filtered_data:
        g = d['seq'].count('G')
        c = d['seq'].count('C')
        tot = len(d['seq'])
        gc = (g + c) / tot * 100 if tot > 0 else 0
        skew = (g - c) / (g + c) if (g + c) > 0 else 0
        table_rows.append({
            'id': d['id'][:15],
            'length': tot,
            'gc': round(gc, 1),
            'at': round(100 - gc, 1),
            'skew': round(skew, 2),
            'entropy': hitung_entropi(d['seq']),
            'orf_len': math.floor(tot / 3)
        })

    return jsonify({
        'summary': {
            'gc_percent': round(gc_percent, 1),
            'start_found': start_found,
            'stop_found': stop_found,
            'total_bp': total_basa
        },
        'viewer': {
            'dna': seq[:80] + "...",
            'mrna': mrna[:80] + "...",
            'peptide': peptide[:80] + "..."
        },
        'charts': {
            'skew_labels': skew_labels,
            'skew_data': gc_skew_data,
            'purines': [a_count + g_count, c_count + t_count],
            'hydro_data': hydro_data,
            'amino_prop': [polar_count, nonpolar_count, start_found, stop_found]
        },
        'table': table_rows
    })

if __name__ == '__main__':
    # Menjalankan server Flask secara lokal
    app.run(debug=True, port=5000)