
from flask import Flask, jsonify, request
from flask_cors import CORS
import requests
import json
import time
import re
from Bio import Entrez
from Bio.Blast import NCBIWWW
from bs4 import BeautifulSoup

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
        
        # Search PDB
        pdb_results = search_pdb(query)
        results.extend(pdb_results)
        
        # Search PubMed
        pubmed_results = search_pubmed(query)
        results.extend(pubmed_results)
        
        # Search UCSC Genome Browser
        ucsc_results = search_ucsc(query)
        results.extend(ucsc_results)
        
        # Search DDBJ
        ddbj_results = search_ddbj(query)
        results.extend(ddbj_results)
        
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
    # First, check if it's in FASTA format
    if query.startswith('>'):
        lines = query.strip().split('\n')
        if len(lines) > 1:
            # Extract just the sequence part
            sequence = ''.join(lines[1:])
            query = sequence
    
    # Check if the sequence consists primarily of valid nucleotide or amino acid characters
    nucleotide_chars = set('ATCGNU atcgnu')
    amino_acid_chars = set('ACDEFGHIKLMNPQRSTVWY acdefghiklmnpqrstvwy')
    
    # Count valid characters
    nucleotide_count = sum(1 for c in query if c in nucleotide_chars)
    amino_acid_count = sum(1 for c in query if c in amino_acid_chars)
    
    # If at least 80% of the characters are valid nucleotides/amino acids
    return (nucleotide_count / len(query) > 0.8 or amino_acid_count / len(query) > 0.8) and len(query) >= 10

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

def search_pdb(query):
    try:
        # First check if it's a direct PDB ID (4 characters, alphanumeric)
        if re.match(r'^[a-zA-Z0-9]{4}$', query):
            return [{
                "id": f"pdb-{query}",
                "title": f"PDB Structure: {query}",
                "description": "Protein Data Bank structure",
                "source": "pdb",
                "url": f"https://www.rcsb.org/structure/{query}",
                "additionalData": {
                    "id": query,
                    "type": "Structure"
                }
            }]
        
        # Otherwise search the PDB
        response = requests.get(
            f"https://search.rcsb.org/rcsbsearch/v2/query?json={{\"query\":{{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{{\"value\":\"{query}\"}}}}}}",
        )
        
        if response.ok:
            data = response.json()
            results = []
            
            for item in data.get('result_set', [])[:5]:
                pdb_id = item.get('identifier', '')
                
                # Get additional data
                detail_response = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
                detail_data = {}
                if detail_response.ok:
                    detail_data = detail_response.json()
                
                results.append({
                    "id": f"pdb-{pdb_id}",
                    "title": f"PDB: {detail_data.get('struct', {}).get('title', pdb_id)}",
                    "description": detail_data.get('struct', {}).get('pdbx_descriptor', 'Protein structure'),
                    "source": "pdb",
                    "url": f"https://www.rcsb.org/structure/{pdb_id}",
                    "additionalData": {
                        "resolution": detail_data.get('rcsb_entry_info', {}).get('resolution_combined', 'Unknown'),
                        "experimental_method": detail_data.get('exptl', {}).get('method', 'Unknown'),
                        "release_date": detail_data.get('rcsb_accession_info', {}).get('deposit_date', 'Unknown')
                    }
                })
            
            return results
        return []
    except Exception as e:
        print(f"PDB search error: {e}")
        return []

def search_pubmed(query):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
        record = Entrez.read(handle)
        
        results = []
        for pubmed_id in record.get("IdList", []):
            # Fetch article details
            article_handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
            article_data = Entrez.read(article_handle)
            
            if article_data.get("PubmedArticle", []):
                article = article_data["PubmedArticle"][0]
                article_details = article.get("MedlineCitation", {}).get("Article", {})
                
                # Extract title and abstract
                title = article_details.get("ArticleTitle", "No title available")
                
                abstract_text = "No abstract available"
                if "Abstract" in article_details and "AbstractText" in article_details["Abstract"]:
                    abstract_parts = article_details["Abstract"]["AbstractText"]
                    if isinstance(abstract_parts, list):
                        abstract_text = " ".join([str(part) for part in abstract_parts])
                    else:
                        abstract_text = str(abstract_parts)
                
                # Extract journal and publication date
                journal = article_details.get("Journal", {}).get("Title", "Unknown Journal")
                
                pub_date = "Unknown Date"
                if "PubDate" in article_details.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {}):
                    pub_date_parts = article_details["Journal"]["JournalIssue"]["PubDate"]
                    pub_date = " ".join([pub_date_parts.get(k, "") for k in ["Year", "Month", "Day"] if k in pub_date_parts])
                
                results.append({
                    "id": f"pubmed-{pubmed_id}",
                    "title": title,
                    "description": abstract_text[:250] + "..." if len(abstract_text) > 250 else abstract_text,
                    "source": "pubmed",
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/",
                    "additionalData": {
                        "journal": journal,
                        "publication_date": pub_date,
                        "pmid": pubmed_id
                    }
                })
            
        return results
    except Exception as e:
        print(f"PubMed search error: {e}")
        return []

def search_ucsc(query):
    try:
        # For UCSC, we'll check if it's a gene or genome region
        if re.match(r'^chr[0-9XY]+:[0-9]+-[0-9]+$', query):  # Matches patterns like chr1:12345-67890
            # It's a genome region
            chrom, pos = query.split(':')
            start, end = pos.split('-')
            
            return [{
                "id": f"ucsc-{query}",
                "title": f"Genome Region: {query}",
                "description": f"Genomic region coordinates: {chrom}:{start}-{end}",
                "source": "ucsc",
                "url": f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={query}",
                "additionalData": {
                    "chromosome": chrom,
                    "start": start,
                    "end": end
                }
            }]
        else:
            # Assume it's a gene symbol
            return [{
                "id": f"ucsc-{query}",
                "title": f"UCSC Gene: {query}",
                "description": f"UCSC Genome Browser entry for gene {query}",
                "source": "ucsc",
                "url": f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={query}",
                "additionalData": {
                    "gene": query,
                    "assembly": "hg38"
                }
            }]
    except Exception as e:
        print(f"UCSC search error: {e}")
        return []

def search_ddbj(query):
    try:
        # DDBJ doesn't have a simple REST API, so we'll make a basic request to their search
        # This is a simplified approach - in production you'd want more robust parsing
        url = f"https://getentry.ddbj.nig.ac.jp/getentry/na/{query}?filetype=html"
        response = requests.get(url)
        
        results = []
        
        if response.ok and 'No such entry' not in response.text:
            # Try to parse some basic information
            soup = BeautifulSoup(response.text, 'lxml')
            title_elem = soup.find('title')
            title = title_elem.text if title_elem else f"DDBJ: {query}"
            
            # Get description if available
            description = "Nucleotide sequence data from DDBJ"
            definition_elem = soup.find('span', class_='definition')
            if definition_elem:
                description = definition_elem.text
            
            results.append({
                "id": f"ddbj-{query}",
                "title": title,
                "description": description,
                "source": "ddbj",
                "url": f"https://getentry.ddbj.nig.ac.jp/getentry/na/{query}",
                "additionalData": {
                    "accession": query
                }
            })
        
        return results
    except Exception as e:
        print(f"DDBJ search error: {e}")
        return []

def blast_search(query):
    try:
        print("Running BLAST search...")
        # Remove FASTA header if present
        if query.startswith('>'):
            lines = query.strip().split('\n')
            if len(lines) > 1:
                sequence = ''.join(lines[1:])
            else:
                sequence = ''
        else:
            sequence = query
        
        # Determine if it's nucleotide or protein
        nucleotide_chars = set('ATCGNU atcgnu')
        amino_acid_chars = set('ACDEFGHIKLMNPQRSTVWY acdefghiklmnpqrstvwy')
        
        nucleotide_count = sum(1 for c in sequence if c in nucleotide_chars)
        amino_acid_count = sum(1 for c in sequence if c in amino_acid_chars)
        
        if nucleotide_count / len(sequence) > amino_acid_count / len(sequence):
            program = "blastn"
            database = "nt"
        else:
            program = "blastp"
            database = "nr"
        
        # For longer sequences, use a shorter expectation value to speed up the search
        evalue = 1e-3 if len(sequence) > 100 else 10
        
        result_handle = NCBIWWW.qblast(program, database, sequence, expect=evalue)
        
        # Process the results to extract useful information
        blast_results = result_handle.read()
        
        # We'll return this as a special result
        return [{
            "id": "blast-result",
            "title": f"BLAST {program.upper()} Results",
            "description": f"Sequence alignment results against {database} database.",
            "source": "ncbi",
            "url": "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
            "additionalData": {
                "program": program,
                "database": database,
                "result": blast_results[:1000] + "..." if len(blast_results) > 1000 else blast_results
            }
        }]
    except Exception as e:
        print(f"BLAST search error: {e}")
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
        
        # Check if we have PDB results (protein structure)
        if 'pdb' in sources:
            count = len(sources['pdb'])
            result = sources['pdb'][0]
            summary_parts.append(f"Found {count} protein structure(s) in PDB.")
            if 'additionalData' in result and 'resolution' in result['additionalData']:
                summary_parts.append(f"Best resolution: {result['additionalData']['resolution']}.")
        
        # Check for gene information
        if 'ensembl' in sources:
            result = sources['ensembl'][0]
            summary_parts.append(f"Found gene {result['title']} in Ensembl database.")
            if 'additionalData' in result:
                if 'chromosome' in result['additionalData']:
                    summary_parts.append(f"Located on chromosome {result['additionalData']['chromosome']}.")
                if 'type' in result['additionalData']:
                    summary_parts.append(f"Gene type: {result['additionalData']['type']}.")
        
        if 'ncbi' in sources:
            result = sources['ncbi'][0]
            summary_parts.append(f"NCBI Gene Database entry: {result['title']}.")
            if 'additionalData' in result and 'organism' in result['additionalData']:
                summary_parts.append(f"Organism: {result['additionalData']['organism']}.")
        
        # Check for protein information
        if 'uniprot' in sources:
            result = sources['uniprot'][0]
            summary_parts.append(f"Associated protein in UniProt: {result['title']}.")
            if 'additionalData' in result and 'length' in result['additionalData']:
                summary_parts.append(f"Protein length: {result['additionalData']['length']} amino acids.")
        
        # Check for sequence data
        if 'ebi' in sources:
            count = len(sources['ebi'])
            summary_parts.append(f"Found {count} related sequence(s) in EBI ENA.")
        
        if 'ddbj' in sources:
            count = len(sources['ddbj'])
            summary_parts.append(f"Found {count} nucleotide sequence(s) in DDBJ database.")
        
        # Check for literature
        if 'pubmed' in sources:
            count = len(sources['pubmed'])
            summary_parts.append(f"Found {count} relevant publication(s) in PubMed.")
            if count > 0:
                result = sources['pubmed'][0]
                if 'additionalData' in result and 'journal' in result['additionalData']:
                    summary_parts.append(f"Most recent in {result['additionalData']['journal']}.")
        
        # Check for GitHub repositories
        if 'github' in sources:
            count = len(sources['github'])
            summary_parts.append(f"Found {count} related bioinformatics repositories on GitHub.")
            if count > 0:
                stars = sources['github'][0]['additionalData'].get('stars', 0)
                summary_parts.append(f"Most popular has {stars} stars.")
        
        # Check for BLAST results
        if any(result['id'] == 'blast-result' for result in results):
            summary_parts.append("BLAST search completed. View alignment results for details.")
        
        summary = " ".join(summary_parts)
        if not summary:
            summary = "Found results across multiple biological databases. Check individual sections for details."
        
        return jsonify({"summary": summary})
    except Exception as e:
        print(f"Summary generation error: {e}")
        return jsonify({"summary": "An error occurred while generating the summary."})

if __name__ == '__main__':
    app.run(debug=True, port=5000)
