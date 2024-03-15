import io
import json
from Bio import SeqIO
from flask import Flask, render_template_string, request, jsonify

app = Flask(__name__)

def find_conserved_regions(proteins, reference_seq):
    conserved_regions = []
    for i, protein in enumerate(proteins):
        if i == 0:
            continue
        protein_len = len(protein)
        reference_len = len(reference_seq)
        start = 0
        while start < min(protein_len, reference_len):
            if protein[start] == reference_seq[start]:
                end = start + 1
                while end < min(protein_len, reference_len) and protein[end] == reference_seq[end]:
                    end += 1
                conserved_regions.append((i, start, end))
                start = end
            else:
                start += 1
    return conserved_regions

@app.route("/")
def index():
    return render_template_string("""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Conserved Regions</title>
    </head>
    <body>
        <h1>Upload a FASTA file</h1>
        <form action="{{ url_for('run_script') }}" method="post" enctype="multipart/form-data">
            <input type="file" name="file" accept=".fasta" required>
            <button type="submit">Run</button>
        </form>
    </body>
    </html>
    """)

@app.route('/run-script', methods=['POST'])
def run_script():
    if request.method == 'POST':
        file_content = request.files['file'].read().decode('utf-8')
        if file_content:
            proteins = []
            for record in SeqIO.parse(io.StringIO(file_content), 'fasta'):
                proteins.append(str(record.seq))

            reference_seq = proteins[0]
            del proteins[0]

            conserved_regions = find_conserved_regions(proteins, reference_seq)

            output = {
                "number_of_proteins": len(proteins),
                "number_of_conserved_regions": len(conserved_regions),
                "conserved_regions": []
            }

            for i, protein in enumerate(proteins):
                output["conserved_regions"].append({
                    "protein_index": i + 1,
                    "conserved_regions": conserved_regions[i]
                })

            return jsonify(output)

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8080, debug=True)