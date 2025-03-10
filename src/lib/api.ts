import { ResultItem, DatabaseSource } from "../components/ResultCard";

// Simulate API delay for development/demo
const simulateApiDelay = () => new Promise(resolve => setTimeout(resolve, Math.random() * 1000 + 500));

// Example data for mock results
const mockResults: Record<string, ResultItem[]> = {
  "BRCA1": [
    {
      id: "ENSG00000012048",
      title: "BRCA1 (ENSG00000012048)",
      description: "DNA repair associated. Plays a role in maintaining genomic stability and acts as a tumor suppressor.",
      source: "ensembl",
      url: "https://ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000012048",
      additionalData: {
        species: "Homo sapiens",
        chromosome: "17",
        type: "Protein coding"
      }
    },
    {
      id: "672",
      title: "BRCA1 - breast cancer type 1 susceptibility protein",
      description: "This gene provides instructions for making a protein that acts as a tumor suppressor.",
      source: "ncbi",
      url: "https://www.ncbi.nlm.nih.gov/gene/672",
      additionalData: {
        organism: "Homo sapiens",
        location: "17q21.31"
      }
    },
    {
      id: "P38398",
      title: "Breast cancer type 1 susceptibility protein",
      description: "E3 ubiquitin-protein ligase that specifically mediates the formation of 'Lys-6'-linked polyubiquitin chains.",
      source: "uniprot",
      url: "https://www.uniprot.org/uniprot/P38398",
      additionalData: {
        length: "1863 aa",
        mass: "207.7 kDa"
      }
    },
    {
      id: "ENA-KU892149",
      title: "Homo sapiens BRCA1, DNA repair associated",
      description: "Complete BRCA1 gene, including promoter region and all exons and introns.",
      source: "ebi",
      url: "https://www.ebi.ac.uk/ena/browser/view/KU892149",
      additionalData: {
        length: "81189 bp"
      }
    },
    {
      id: "DDBJ-AB709907",
      title: "Human BRCA1 gene sequence",
      description: "Homo sapiens BRCA1 mRNA, partial sequence, alternative transcript.",
      source: "ddbj",
      url: "https://www.ddbj.nig.ac.jp/getentry?database=na&accession_number=AB709907",
      additionalData: {
        released: "2013-03-15"
      }
    },
    {
      id: "github-brca1-analysis",
      title: "BRCA1-mutation-analysis",
      description: "Repository containing tools for analyzing BRCA1 mutations and their effects on protein structure.",
      source: "github",
      url: "https://github.com/example/BRCA1-mutation-analysis",
      additionalData: {
        stars: "124",
        language: "Python"
      }
    }
  ],
  "p53": [
    {
      id: "ENSG00000141510",
      title: "TP53 (ENSG00000141510)",
      description: "Tumor protein p53, acts as a tumor suppressor in many tumor types; induces growth arrest or apoptosis.",
      source: "ensembl",
      url: "https://ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000141510",
      additionalData: {
        species: "Homo sapiens",
        chromosome: "17",
        type: "Protein coding"
      }
    },
    {
      id: "7157",
      title: "TP53 - tumor protein p53",
      description: "This gene encodes a tumor suppressor protein containing transcriptional activation, DNA binding, and oligomerization domains.",
      source: "ncbi",
      url: "https://www.ncbi.nlm.nih.gov/gene/7157",
      additionalData: {
        organism: "Homo sapiens",
        location: "17p13.1"
      }
    },
    {
      id: "P04637",
      title: "Cellular tumor antigen p53",
      description: "Acts as a tumor suppressor in many tumor types; induces growth arrest or apoptosis depending on the physiological circumstances and cell type.",
      source: "uniprot",
      url: "https://www.uniprot.org/uniprot/P04637",
      additionalData: {
        length: "393 aa",
        mass: "43.7 kDa"
      }
    }
  ],
  "insulin": [
    {
      id: "ENSG00000254647",
      title: "INS (ENSG00000254647)",
      description: "Insulin, produced in the pancreas, is central to regulating carbohydrate and fat metabolism in the body.",
      source: "ensembl",
      url: "https://ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000254647",
      additionalData: {
        species: "Homo sapiens",
        chromosome: "11",
        type: "Protein coding"
      }
    },
    {
      id: "3630",
      title: "INS - insulin",
      description: "This gene encodes insulin, a peptide hormone that regulates carbohydrate and fat metabolism in the body.",
      source: "ncbi",
      url: "https://www.ncbi.nlm.nih.gov/gene/3630",
      additionalData: {
        organism: "Homo sapiens",
        location: "11p15.5"
      }
    }
  ],
  "ATGGCGCGAT": [
    {
      id: "ENA-MK341985",
      title: "Synthetic construct sequence containing ATGGCGCGAT",
      description: "Synthetic DNA construct containing multiple restriction sites and cloning sequences.",
      source: "ebi",
      url: "https://www.ebi.ac.uk/ena/browser/view/MK341985",
      additionalData: {
        length: "3025 bp"
      }
    },
    {
      id: "DDBJ-LC123456",
      title: "DNA sequence with ATGGCGCGAT motif",
      description: "Partial sequence containing conserved ATGGCGCGAT motif found in multiple species.",
      source: "ddbj",
      url: "https://www.ddbj.nig.ac.jp/getentry?database=na&accession_number=LC123456",
      additionalData: {
        released: "2022-01-10"
      }
    }
  ],
  "coronavirus": [
    {
      id: "ENSG00000271049",
      title: "SARS-CoV-2 related sequence in human",
      description: "Genetic element related to coronavirus integration in human genome.",
      source: "ensembl",
      url: "https://ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000271049",
      additionalData: {
        species: "Homo sapiens",
        type: "Non-coding"
      }
    },
    {
      id: "1489680",
      title: "SARS-CoV-2 (COVID-19 virus)",
      description: "Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome.",
      source: "ncbi",
      url: "https://www.ncbi.nlm.nih.gov/nuccore/NC_045512",
      additionalData: {
        length: "29903 bp",
        type: "Reference genome"
      }
    },
    {
      id: "github-covid-analysis",
      title: "coronavirus-genome-analysis",
      description: "Tools for analyzing coronavirus genomes, mutations, and evolution patterns.",
      source: "github",
      url: "https://github.com/example/coronavirus-genome-analysis",
      additionalData: {
        stars: "2135",
        language: "R"
      }
    }
  ]
};

// Default fallback search results
const defaultResults: ResultItem[] = [
  {
    id: "default-1",
    title: "Unknown query result",
    description: "Your search didn't match any specific examples in our mock database. In a real implementation, we would query actual bioinformatics databases.",
    source: "ncbi",
    url: "https://www.ncbi.nlm.nih.gov/",
    additionalData: {
      note: "Mock result"
    }
  },
  {
    id: "default-2",
    title: "Try a sample query",
    description: "Try searching for examples like 'BRCA1', 'p53', 'insulin', 'ATGGCGCGAT', or 'coronavirus'.",
    source: "ensembl",
    url: "https://ensembl.org/",
    additionalData: {
      examples: "BRCA1, p53, insulin"
    }
  }
];

// Mock search function that will be replaced with real API calls in production
export async function searchBioDatabases(query: string): Promise<ResultItem[]> {
  // For development purposes, we'll use mock data with a simulated delay
  // In production, this would call the actual APIs below
  await simulateApiDelay();
  
  // For demo purposes, we'll match some predefined queries
  const normalizedQuery = query.trim().toLowerCase();
  
  // Check for exact matches first
  if (mockResults[query]) {
    return mockResults[query];
  }
  
  // Check for partial matches
  for (const [key, results] of Object.entries(mockResults)) {
    if (key.toLowerCase().includes(normalizedQuery) || normalizedQuery.includes(key.toLowerCase())) {
      return results;
    }
  }
  
  // Return default results if no match
  return defaultResults;
}

// Example implementation to search Ensembl for gene data
async function searchEnsembl(term: string): Promise<ResultItem[]> {
  try {
    const response = await fetch(`https://rest.ensembl.org/lookup/symbol/homo_sapiens/${term}?content-type=application/json`);
    const data = await response.json();
    
    if (data.id) {
      return [{
        id: data.id,
        title: `${data.display_name} (${data.id})`,
        description: data.description || "No description available",
        source: "ensembl",
        url: `https://ensembl.org/Homo_sapiens/Gene/Summary?g=${data.id}`,
        additionalData: {
          species: data.species || "Homo sapiens",
          chromosome: data.seq_region_name || "Unknown",
          type: data.biotype || "Unknown"
        }
      }];
    }
    return [];
  } catch (error) {
    console.error("Ensembl API error:", error);
    return [];
  }
}

// Example implementation to search NCBI for gene data
async function searchNCBI(term: string): Promise<ResultItem[]> {
  try {
    const response = await fetch(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${term}&retmode=json&retmax=5`);
    const data = await response.json();
    
    if (data.esearchresult && data.esearchresult.idlist && data.esearchresult.idlist.length > 0) {
      return data.esearchresult.idlist.map((id: string) => ({
        id,
        title: `${term} (${id})`,
        description: "NCBI Gene Database Entry",
        source: "ncbi" as DatabaseSource,
        url: `https://www.ncbi.nlm.nih.gov/gene/${id}`,
        additionalData: {
          organism: "Homo sapiens", // Would be fetched from a summary API in production
          id: id
        }
      }));
    }
    return [];
  } catch (error) {
    console.error("NCBI API error:", error);
    return [];
  }
}

// Example implementation to search UniProt for protein data
async function searchUniProt(term: string): Promise<ResultItem[]> {
  try {
    const response = await fetch(`https://rest.uniprot.org/uniprotkb/search?query=${term}&size=5&format=json`);
    const data = await response.json();
    
    if (data.results && data.results.length > 0) {
      return data.results.map((result: any) => {
        const name = result.proteinDescription?.recommendedName?.fullName?.value || 'No name';
        return {
          id: result.primaryAccession,
          title: name,
          description: result.comments?.find((c: any) => c.commentType === 'FUNCTION')?.texts[0]?.value || "No description available",
          source: "uniprot" as DatabaseSource,
          url: `https://www.uniprot.org/uniprotkb/${result.primaryAccession}`,
          additionalData: {
            length: result.sequence?.length ? `${result.sequence.length} aa` : "Unknown",
            mass: result.sequence?.molWeight ? `${(result.sequence.molWeight/1000).toFixed(1)} kDa` : "Unknown"
          }
        };
      });
    }
    return [];
  } catch (error) {
    console.error("UniProt API error:", error);
    return [];
  }
}

// Example implementation to search EBI ENA for sequence data
async function searchENA(term: string): Promise<ResultItem[]> {
  try {
    const response = await fetch(`https://www.ebi.ac.uk/ena/portal/api/search?query=${term}&result=sequence&limit=5`);
    const data = await response.json();
    
    if (data && data.length > 0) {
      return data.map((hit: any) => ({
        id: hit.accession,
        title: `${hit.description || 'ENA Sequence'} (${hit.accession})`,
        description: hit.description || "No description available",
        source: "ebi" as DatabaseSource,
        url: `https://www.ebi.ac.uk/ena/browser/view/${hit.accession}`,
        additionalData: {
          length: hit.length ? `${hit.length} bp` : "Unknown"
        }
      }));
    }
    return [];
  } catch (error) {
    console.error("EBI ENA API error:", error);
    return [];
  }
}

// Example implementation to search DDBJ for sequence data
async function searchDDBJ(term: string): Promise<ResultItem[]> {
  // Note: This is a placeholder. DDBJ may not have a public API like this.
  // A server-side proxy might be needed for production.
  try {
    // In a real implementation, this would use the actual DDBJ API endpoint
    // const response = await fetch(`https://ddbj.nig.ac.jp/api/search?query=${term}`);
    // const data = await response.json();
    
    // For now, returning mock DDBJ data based on term
    if (mockResults[term]?.filter(r => r.source === "ddbj").length > 0) {
      return mockResults[term].filter(r => r.source === "ddbj");
    }
    return [];
  } catch (error) {
    console.error("DDBJ API error:", error);
    return [];
  }
}

// Example implementation to search GitHub for bioinformatics repositories
async function searchGitHub(term: string): Promise<ResultItem[]> {
  try {
    const response = await fetch(`https://api.github.com/search/repositories?q=${term}+bioinformatics`);
    const data = await response.json();
    
    if (data.items && data.items.length > 0) {
      return data.items.slice(0, 5).map((item: any) => ({
        id: item.id.toString(),
        title: item.full_name,
        description: item.description || "No description available",
        source: "github" as DatabaseSource,
        url: item.html_url,
        additionalData: {
          stars: item.stargazers_count.toString(),
          language: item.language || "Not specified"
        }
      }));
    }
    return [];
  } catch (error) {
    console.error("GitHub API error:", error);
    return [];
  }
}

// Production-ready search function (would replace mockSearchBioDatabases)
// This is commented out as we're using mock data for the demo
/*
export async function searchBioDatabases(query: string): Promise<ResultItem[]> {
  try {
    // Execute all searches in parallel
    const [ensemblResults, ncbiResults, uniprotResults, enaResults, ddbjResults, githubResults] = await Promise.all([
      searchEnsembl(query),
      searchNCBI(query),
      searchUniProt(query),
      searchENA(query),
      searchDDBJ(query),
      searchGitHub(query)
    ]);
    
    // Combine all results
    return [
      ...ensemblResults,
      ...ncbiResults,
      ...uniprotResults,
      ...enaResults,
      ...ddbjResults,
      ...githubResults
    ];
  } catch (error) {
    console.error("Search error:", error);
    return defaultResults;
  }
}
*/

// Generate summary based on results
export function generateSummary(results: ResultItem[]): string {
  if (results.length === 0) {
    return "No results found. Try refining your search.";
  }
  
  // Prioritize results from certain databases
  const ensemblResult = results.find(r => r.source === 'ensembl');
  const uniprotResult = results.find(r => r.source === 'uniprot');
  const ncbiResult = results.find(r => r.source === 'ncbi');
  
  if (ensemblResult) {
    const additionalInfo = ensemblResult.additionalData;
    return `${ensemblResult.title} - ${additionalInfo?.species || 'Species unknown'}${additionalInfo?.chromosome ? `, Chromosome ${additionalInfo.chromosome}` : ''}. ${ensemblResult.description.split('.')[0]}.`;
  } else if (uniprotResult) {
    return `${uniprotResult.title} - ${uniprotResult.additionalData?.length || ''} ${uniprotResult.additionalData?.mass ? `(${uniprotResult.additionalData.mass})` : ''}. ${uniprotResult.description.split('.')[0]}.`;
  } else if (ncbiResult) {
    return `${ncbiResult.title} - ${ncbiResult.additionalData?.organism || ''}. ${ncbiResult.description.split('.')[0]}.`;
  }
  
  // Fallback to the first result
  return `${results[0].title} - ${results[0].description.split('.')[0]}.`;
}

// Group results by database source
export function groupResultsBySource(results: ResultItem[]): Record<DatabaseSource, ResultItem[]> {
  const grouped: Record<DatabaseSource, ResultItem[]> = {
    ensembl: [],
    ncbi: [],
    uniprot: [],
    ebi: [],
    ddbj: [],
    github: []
  };
  
  results.forEach(result => {
    grouped[result.source].push(result);
  });
  
  return grouped;
}
