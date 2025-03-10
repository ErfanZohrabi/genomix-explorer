
from flask import Flask, jsonify, request
from flask_cors import CORS
import requests
import json

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Default fallback search results
default_results = [
    {
        "id": "default-1",
        "title": "Unknown query result",
        "description": "Your search didn't match any specific examples in our mock database. In a real implementation, we would query actual bioinformatics databases.",
        "source": "ncbi",
        "url": "https://www.ncbi.nlm.nih.gov/",
        "additionalData": {
            "note": "Mock result"
        }
    },
    {
        "id": "default-2",
        "title": "Try a sample query",
        "description": "Try searching for examples like 'BRCA1', 'p53', 'insulin', 'ATGGCGCGAT', or 'coronavirus'.",
        "source": "ensembl",
        "url": "https://ensembl.org/",
        "additionalData": {
            "examples": "BRCA1, p53, insulin"
        }
    }
]

# Search Ensembl for gene data
def search_ensembl(term):
    try:
        response = requests.get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{term}?content-type=application/json")
        data = response.json()
        
        if 'id' in data:
            return [{
                "id": data['id'],
                "title": f"{data.get('display_name', term)} ({data['id']})",
                "description": data.get('description', "No description available"),
                "source": "ensembl",
                "url": f"https://ensembl.org/Homo_sapiens/Gene/Summary?g={data['id']}",
                "additionalData": {
                    "species": data.get('species', "Homo sapiens"),
                    "chromosome": data.get('seq_region_name', "Unknown"),
                    "type": data.get('biotype', "Unknown")
                }
            }]
        return []
    except Exception as e:
        print(f"Ensembl API error: {e}")
        return []

# Search NCBI for gene data
def search_ncbi(term):
    try:
        response = requests.get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={term}&retmode=json&retmax=5")
        data = response.json()
        
        if 'esearchresult' in data and 'idlist' in data['esearchresult'] and data['esearchresult']['idlist']:
            return [{
                "id": gene_id,
                "title": f"{term} ({gene_id})",
                "description": "NCBI Gene Database Entry",
                "source": "ncbi",
                "url": f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}",
                "additionalData": {
                    "organism": "Homo sapiens",  # Would be fetched from a summary API in production
                    "id": gene_id
                }
            } for gene_id in data['esearchresult']['idlist']]
        return []
    except Exception as e:
        print(f"NCBI API error: {e}")
        return []

# Search UniProt for protein data
def search_uniprot(term):
    try:
        response = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={term}&size=5&format=json")
        data = response.json()
        
        if 'results' in data and data['results']:
            results = []
            for result in data['results']:
                name = "No name"
                if 'proteinDescription' in result and 'recommendedName' in result['proteinDescription']:
                    if 'fullName' in result['proteinDescription']['recommendedName']:
                        name = result['proteinDescription']['recommendedName']['fullName'].get('value', 'No name')
                
                description = "No description available"
                if 'comments' in result:
                    for comment in result['comments']:
                        if comment['commentType'] == 'FUNCTION' and 'texts' in comment and comment['texts']:
                            description = comment['texts'][0].get('value', description)
                
                results.append({
                    "id": result['primaryAccession'],
                    "title": name,
                    "description": description,
                    "source": "uniprot",
                    "url": f"https://www.uniprot.org/uniprotkb/{result['primaryAccession']}",
                    "additionalData": {
                        "length": f"{result['sequence']['length']} aa" if 'sequence' in result and 'length' in result['sequence'] else "Unknown",
                        "mass": f"{(result['sequence']['molWeight']/1000):.1f} kDa" if 'sequence' in result and 'molWeight' in result['sequence'] else "Unknown"
                    }
                })
            return results
        return []
    except Exception as e:
        print(f"UniProt API error: {e}")
        return []

# Search EBI ENA for sequence data
def search_ena(term):
    try:
        response = requests.get(f"https://www.ebi.ac.uk/ena/portal/api/search?query={term}&result=sequence&limit=5")
        data = response.json()
        
        if data:
            return [{
                "id": hit['accession'],
                "title": f"{hit.get('description', 'ENA Sequence')} ({hit['accession']})",
                "description": hit.get('description', "No description available"),
                "source": "ebi",
                "url": f"https://www.ebi.ac.uk/ena/browser/view/{hit['accession']}",
                "additionalData": {
                    "length": f"{hit['length']} bp" if 'length' in hit else "Unknown"
                }
            } for hit in data]
        return []
    except Exception as e:
        print(f"EBI ENA API error: {e}")
        return []

# Search DDBJ for sequence data (placeholder - requires implementation based on actual API)
def search_ddbj(term):
    # Note: This is a placeholder. DDBJ may not have a public API like this.
    # A server-side proxy might be needed for production.
    try:
        # In a real implementation, this would use the actual DDBJ API endpoint
        # For now, returning an empty list
        return []
    except Exception as e:
        print(f"DDBJ API error: {e}")
        return []

# Search GitHub for bioinformatics repositories
def search_github(term):
    try:
        response = requests.get(f"https://api.github.com/search/repositories?q={term}+bioinformatics")
        data = response.json()
        
        if 'items' in data and data['items']:
            return [{
                "id": str(item['id']),
                "title": item['full_name'],
                "description": item.get('description', "No description available"),
                "source": "github",
                "url": item['html_url'],
                "additionalData": {
                    "stars": str(item['stargazers_count']),
                    "language": item.get('language', "Not specified")
                }
            } for item in data['items'][:5]]  # Limit to first 5 results
        return []
    except Exception as e:
        print(f"GitHub API error: {e}")
        return []

@app.route('/api/search', methods=['GET'])
def search_databases():
    term = request.args.get('query', '')
    if not term:
        return jsonify(default_results)
    
    # Run all searches in parallel (for a production app, you would use async)
    ensembl_results = search_ensembl(term)
    ncbi_results = search_ncbi(term)
    uniprot_results = search_uniprot(term)
    ena_results = search_ena(term)
    ddbj_results = search_ddbj(term)
    github_results = search_github(term)
    
    # Combine all results
    all_results = ensembl_results + ncbi_results + uniprot_results + ena_results + ddbj_results + github_results
    
    # Return default results if no matches found
    if not all_results:
        return jsonify(default_results)
    
    return jsonify(all_results)

@app.route('/api/summary', methods=['POST'])
def generate_summary():
    results = request.json.get('results', [])
    
    if not results:
        return jsonify({"summary": "No results found. Try refining your search."})
    
    # Prioritize results from certain databases
    ensembl_result = next((r for r in results if r['source'] == 'ensembl'), None)
    uniprot_result = next((r for r in results if r['source'] == 'uniprot'), None)
    ncbi_result = next((r for r in results if r['source'] == 'ncbi'), None)
    
    if ensembl_result:
        additional_info = ensembl_result.get('additionalData', {})
        species = additional_info.get('species', 'Species unknown')
        chromosome = additional_info.get('chromosome', '')
        chrom_text = f", Chromosome {chromosome}" if chromosome else ''
        description = ensembl_result.get('description', '').split('.')[0]
        return jsonify({"summary": f"{ensembl_result['title']} - {species}{chrom_text}. {description}."})
    
    if uniprot_result:
        additional_info = uniprot_result.get('additionalData', {})
        length = additional_info.get('length', '')
        mass = additional_info.get('mass', '')
        mass_text = f" ({mass})" if mass else ''
        description = uniprot_result.get('description', '').split('.')[0]
        return jsonify({"summary": f"{uniprot_result['title']} - {length}{mass_text}. {description}."})
    
    if ncbi_result:
        additional_info = ncbi_result.get('additionalData', {})
        organism = additional_info.get('organism', '')
        description = ncbi_result.get('description', '').split('.')[0]
        return jsonify({"summary": f"{ncbi_result['title']} - {organism}. {description}."})
    
    # Fallback to the first result
    description = results[0].get('description', '').split('.')[0]
    return jsonify({"summary": f"{results[0]['title']} - {description}."})

if __name__ == '__main__':
    app.run(debug=True, port=5000)
