
from flask import Flask, jsonify, request
from flask_cors import CORS
import requests
import json
import time
from Bio import Entrez
from Bio.Blast import NCBIWWW

app = Flask(__name__)
CORS(app)

# Configure Entrez
Entrez.email = "your-email@domain.com"  # Replace with your email

@app.route('/api/search', methods=['GET'])
def search_databases():
    query = request.args.get('query', '')
    if not query:
        return jsonify([])
    
    try:
        results = []
        
        # Search NCBI
        ncbi_results = search_ncbi(query)
        results.extend(ncbi_results)
        
        # Search Ensembl
        ensembl_results = search_ensembl(query)
        results.extend(ensembl_results)
        
        # Search UniProt
        uniprot_results = search_uniprot(query)
        results.extend(uniprot_results)
        
        # Search EBI ENA
        ena_results = search_ena(query)
        results.extend(ena_results)
        
        # Search GitHub
        github_results = search_github(query)
        results.extend(github_results)
        
        # Add BLAST search if query looks like a sequence
        if is_sequence(query):
            blast_results = blast_search(query)
            results.extend(blast_results)
        
        return jsonify(results)
    except Exception as e:
        print(f"Search error: {e}")
        return jsonify([])

def is_sequence(query):
    """Check if query looks like a DNA/protein sequence"""
    valid_chars = set('ATCGNU atcgnu')
    return all(c in valid_chars for c in query) and len(query) >= 10

def search_ncbi(query):
    try:
        # Search gene database
        handle = Entrez.esearch(db="gene", term=query, retmax=5)
        record = Entrez.read(handle)
        
        results = []
        for gene_id in record["IdList"]:
            gene_handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="xml")
            gene_record = Entrez.read(gene_handle)
            
            if gene_record:
                results.append({
                    "id": f"ncbi-{gene_id}",
                    "title": f"Gene: {query}",
                    "description": str(gene_record[0].get('description', 'No description available')),
                    "source": "ncbi",
                    "url": f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}",
                    "additionalData": {
                        "organism": gene_record[0].get('organism', {}).get('scientificName', 'Unknown'),
                        "location": gene_record[0].get('locationHist', 'Unknown')
                    }
                })
        
        return results
    except Exception as e:
        print(f"NCBI search error: {e}")
        return []

def search_ensembl(query):
    try:
        response = requests.get(
            f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{query}",
            headers={"Content-Type": "application/json"}
        )
        
        if response.ok:
            data = response.json()
            return [{
                "id": f"ensembl-{data['id']}",
                "title": data.get('display_name', query),
                "description": data.get('description', 'No description available'),
                "source": "ensembl",
                "url": f"https://ensembl.org/Homo_sapiens/Gene/Summary?g={data['id']}",
                "additionalData": {
                    "species": "Homo sapiens",
                    "chromosome": data.get('seq_region_name', 'Unknown'),
                    "type": data.get('biotype', 'Unknown')
                }
            }]
        return []
    except Exception as e:
        print(f"Ensembl search error: {e}")
        return []

def search_uniprot(query):
    try:
        response = requests.get(
            f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=json"
        )
        
        if response.ok:
            data = response.json()
            results = []
            
            for item in data.get('results', [])[:5]:
                name = "Unknown"
                if 'proteinDescription' in item:
                    if 'recommendedName' in item['proteinDescription']:
                        name = item['proteinDescription']['recommendedName'].get('fullName', {}).get('value', 'Unknown')
                
                results.append({
                    "id": f"uniprot-{item['primaryAccession']}",
                    "title": name,
                    "description": get_uniprot_description(item),
                    "source": "uniprot",
                    "url": f"https://www.uniprot.org/uniprotkb/{item['primaryAccession']}",
                    "additionalData": {
                        "accession": item['primaryAccession'],
                        "organism": item.get('organism', {}).get('scientificName', 'Unknown'),
                        "length": item.get('sequence', {}).get('length', 'Unknown')
                    }
                })
            
            return results
        return []
    except Exception as e:
        print(f"UniProt search error: {e}")
        return []

def get_uniprot_description(item):
    if 'comments' in item:
        for comment in item['comments']:
            if comment['commentType'] == 'FUNCTION' and 'texts' in comment:
                return comment['texts'][0].get('value', '')
    return 'No description available'

def search_ena(query):
    try:
        response = requests.get(
            f"https://www.ebi.ac.uk/ena/portal/api/search",
            params={
                "query": query,
                "result": "sequence",
                "limit": 5
            }
        )
        
        if response.ok:
            data = response.json()
            results = []
            
            for item in data:
                results.append({
                    "id": f"ebi-{item['accession']}",
                    "title": item.get('description', 'No title available'),
                    "description": item.get('description', 'No description available'),
                    "source": "ebi",
                    "url": f"https://www.ebi.ac.uk/ena/browser/view/{item['accession']}",
                    "additionalData": {
                        "length": item.get('length', 'Unknown'),
                        "type": item.get('moleculeType', 'Unknown')
                    }
                })
            
            return results
        return []
    except Exception as e:
        print(f"ENA search error: {e}")
        return []

def blast_search(sequence):
    try:
        print("Running BLAST search...")
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        blast_records = result_handle.read()
        
        return [{
            "id": "blast-result",
            "title": "BLAST Search Result",
            "description": "Sequence alignment results",
            "source": "ncbi",
            "url": "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
            "additionalData": {
                "result": str(blast_records)
            }
        }]
    except Exception as e:
        print(f"BLAST search error: {e}")
        return []

def search_github(query):
    try:
        response = requests.get(
            f"https://api.github.com/search/repositories",
            params={
                "q": f"{query} bioinformatics",
                "sort": "stars",
                "order": "desc"
            },
            headers={"Accept": "application/vnd.github.v3+json"}
        )
        
        if response.ok:
            data = response.json()
            results = []
            
            for item in data.get('items', [])[:5]:
                results.append({
                    "id": f"github-{item['id']}",
                    "title": item['full_name'],
                    "description": item.get('description', 'No description available'),
                    "source": "github",
                    "url": item['html_url'],
                    "additionalData": {
                        "stars": item['stargazers_count'],
                        "language": item.get('language', 'Unknown'),
                        "topics": item.get('topics', [])
                    }
                })
            
            return results
        return []
    except Exception as e:
        print(f"GitHub search error: {e}")
        return []

@app.route('/api/summary', methods=['POST'])
def generate_summary():
    try:
        results = request.json.get('results', [])
        
        if not results:
            return jsonify({"summary": "No results found. Try refining your search."})
        
        # Group results by source
        sources = {}
        for result in results:
            source = result['source']
            if source not in sources:
                sources[source] = []
            sources[source].append(result)
        
        # Generate summary
        summary_parts = []
        
        if 'ensembl' in sources:
            result = sources['ensembl'][0]
            summary_parts.append(f"Found gene {result['title']} in Ensembl database.")
        
        if 'uniprot' in sources:
            result = sources['uniprot'][0]
            summary_parts.append(f"Associated protein: {result['title']}.")
        
        if 'ncbi' in sources:
            result = sources['ncbi'][0]
            summary_parts.append(f"NCBI entry: {result['title']}.")
        
        if 'ebi' in sources:
            count = len(sources['ebi'])
            summary_parts.append(f"Found {count} related sequence(s) in EBI ENA.")
        
        if 'github' in sources:
            count = len(sources['github'])
            summary_parts.append(f"Found {count} related bioinformatics repositories on GitHub.")
        
        summary = " ".join(summary_parts)
        if not summary:
            summary = "Found results across multiple biological databases. Check individual sections for details."
        
        return jsonify({"summary": summary})
    except Exception as e:
        print(f"Summary generation error: {e}")
        return jsonify({"summary": "An error occurred while generating the summary."})

if __name__ == '__main__':
    app.run(debug=True, port=5000)

