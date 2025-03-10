
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
    },
    {
      id: "pdb-1jm7",
      title: "BRCA1 RING domain structure",
      description: "Crystal structure of the RING domain from human BRCA1.",
      source: "pdb",
      url: "https://www.rcsb.org/structure/1JM7",
      additionalData: {
        resolution: "2.0 Å",
        experimental_method: "X-ray diffraction",
        release_date: "2001-07-12"
      }
    },
    {
      id: "pubmed-28585547",
      title: "BRCA1 and genome stability",
      description: "This review summarizes the roles of BRCA1 in maintaining genome stability, including DNA repair, cell cycle checkpoint control, and transcriptional regulation.",
      source: "pubmed",
      url: "https://pubmed.ncbi.nlm.nih.gov/28585547/",
      additionalData: {
        journal: "Annual Review of Cancer Biology",
        publication_date: "2017 Mar",
        pmid: "28585547"
      }
    },
    {
      id: "ucsc-chr17:43044295-43125483",
      title: "BRCA1 gene locus",
      description: "Genomic location of the BRCA1 gene on chromosome 17.",
      source: "ucsc",
      url: "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17:43044295-43125483",
      additionalData: {
        chromosome: "chr17",
        start: "43044295",
        end: "43125483"
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
    },
    {
      id: "pdb-1tup",
      title: "Human p53 core domain crystal structure",
      description: "Structure of the core DNA-binding domain of human p53 tumor suppressor protein.",
      source: "pdb",
      url: "https://www.rcsb.org/structure/1tup",
      additionalData: {
        resolution: "2.2 Å",
        experimental_method: "X-ray diffraction",
        release_date: "1994-08-01"
      }
    },
    {
      id: "pubmed-29533785",
      title: "p53: Structure, Function and Therapeutic Applications",
      description: "Review discussing the structure, function, and potential therapeutic applications of p53 in cancer treatment.",
      source: "pubmed",
      url: "https://pubmed.ncbi.nlm.nih.gov/29533785/",
      additionalData: {
        journal: "Journal of Cancer",
        publication_date: "2018 Mar 1",
        pmid: "29533785"
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
    },
    {
      id: "P01308",
      title: "Insulin",
      description: "Insulin stimulates glucose uptake and lipid synthesis and inhibits lipolysis. It is secreted by pancreatic beta cells.",
      source: "uniprot",
      url: "https://www.uniprot.org/uniprot/P01308",
      additionalData: {
        length: "110 aa",
        mass: "11.9 kDa"
      }
    },
    {
      id: "pdb-4ins",
      title: "Human insulin structure",
      description: "Crystal structure of human insulin revealing its hexameric form.",
      source: "pdb",
      url: "https://www.rcsb.org/structure/4ins",
      additionalData: {
        resolution: "1.5 Å",
        experimental_method: "X-ray diffraction",
        release_date: "1989-05-17"
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
    },
    {
      id: "blast-result",
      title: "BLAST Results for ATGGCGCGAT",
      description: "Sequence alignment results for ATGGCGCGAT motif.",
      source: "ncbi",
      url: "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
      additionalData: {
        program: "blastn",
        database: "nt",
        result: "Multiple hits found in various organisms."
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
    },
    {
      id: "pdb-6vxx",
      title: "SARS-CoV-2 spike protein",
      description: "Structure of the SARS-CoV-2 spike glycoprotein in the prefusion conformation.",
      source: "pdb",
      url: "https://www.rcsb.org/structure/6VXX",
      additionalData: {
        resolution: "2.8 Å",
        experimental_method: "Cryo-EM",
        release_date: "2020-03-18"
      }
    },
    {
      id: "pubmed-32015508",
      title: "A Novel Coronavirus from Patients with Pneumonia in China, 2019",
      description: "Initial identification and characterization of the novel coronavirus (later named SARS-CoV-2) from patients with pneumonia in Wuhan, China.",
      source: "pubmed",
      url: "https://pubmed.ncbi.nlm.nih.gov/32015508/",
      additionalData: {
        journal: "New England Journal of Medicine",
        publication_date: "2020 Feb 20",
        pmid: "32015508"
      }
    }
  ],
  "1HHO": [
    {
      id: "pdb-1hho",
      title: "Human Hemoglobin",
      description: "Crystal structure of human deoxyhemoglobin at 1.74 Å resolution.",
      source: "pdb",
      url: "https://www.rcsb.org/structure/1HHO",
      additionalData: {
        resolution: "1.74 Å",
        experimental_method: "X-ray diffraction",
        release_date: "1993-08-11"
      }
    },
    {
      id: "P69905",
      title: "Hemoglobin subunit alpha",
      description: "Part of hemoglobin, responsible for oxygen transport in the blood.",
      source: "uniprot",
      url: "https://www.uniprot.org/uniprot/P69905",
      additionalData: {
        length: "142 aa",
        mass: "15.3 kDa"
      }
    },
    {
      id: "P68871",
      title: "Hemoglobin subunit beta",
      description: "Beta chain of hemoglobin, responsible for oxygen transport.",
      source: "uniprot",
      url: "https://www.uniprot.org/uniprot/P68871",
      additionalData: {
        length: "147 aa",
        mass: "16.0 kDa"
      }
    }
  ],
  "hemoglobin": [
    {
      id: "ENSG00000206172",
      title: "HBA1 (ENSG00000206172)",
      description: "Hemoglobin subunit alpha 1, part of the adult hemoglobin protein.",
      source: "ensembl",
      url: "https://ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000206172",
      additionalData: {
        species: "Homo sapiens",
        chromosome: "16",
        type: "Protein coding"
      }
    },
    {
      id: "P69905",
      title: "Hemoglobin subunit alpha",
      description: "Hemoglobin is the oxygen-carrying protein in red blood cells, composed of two alpha and two beta chains.",
      source: "uniprot",
      url: "https://www.uniprot.org/uniprot/P69905",
      additionalData: {
        length: "142 aa",
        mass: "15.3 kDa"
      }
    },
    {
      id: "pdb-1hho",
      title: "Human Hemoglobin",
      description: "Structure of human deoxyhemoglobin at high resolution.",
      source: "pdb",
      url: "https://www.rcsb.org/structure/1HHO",
      additionalData: {
        resolution: "1.74 Å",
        experimental_method: "X-ray diffraction",
        release_date: "1993-08-11"
      }
    },
    {
      id: "pubmed-31385757",
      title: "Structure and function of hemoglobin",
      description: "Review discussing the structure and function of hemoglobin, including its role in oxygen transport and allosteric regulation.",
      source: "pubmed",
      url: "https://pubmed.ncbi.nlm.nih.gov/31385757/",
      additionalData: {
        journal: "Journal of Biological Chemistry",
        publication_date: "2019 Sep 6",
        pmid: "31385757"
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
    description: "Try searching for examples like 'BRCA1', 'p53', 'insulin', 'ATGGCGCGAT', 'coronavirus', 'hemoglobin', or '1HHO'.",
    source: "ensembl",
    url: "https://ensembl.org/",
    additionalData: {
      examples: "BRCA1, p53, insulin, hemoglobin, 1HHO"
    }
  }
];

// Search function that will call our Python backend API
export async function searchBioDatabases(query: string): Promise<ResultItem[]> {
  try {
    // For development purposes with no backend, use the simulated delay and mock data
    if (process.env.NODE_ENV === 'development' && !process.env.USE_BACKEND) {
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
    
    // Call our Python backend API
    const response = await fetch(`http://localhost:5000/api/search?query=${encodeURIComponent(query)}`);
    if (!response.ok) {
      throw new Error(`API error: ${response.status}`);
    }
    const data = await response.json();
    return data;
  } catch (error) {
    console.error("Search error:", error);
    return defaultResults;
  }
}

// Generate summary based on results using our Python backend
export async function generateSummary(results: ResultItem[]): Promise<string> {
  try {
    // For development purposes with no backend
    if (process.env.NODE_ENV === 'development' && !process.env.USE_BACKEND) {
      if (results.length === 0) {
        return "No results found. Try refining your search.";
      }
      
      // Prioritize results from certain databases
      const ensemblResult = results.find(r => r.source === 'ensembl');
      const uniprotResult = results.find(r => r.source === 'uniprot');
      const ncbiResult = results.find(r => r.source === 'ncbi');
      const pdbResult = results.find(r => r.source === 'pdb');
      const pubmedResults = results.filter(r => r.source === 'pubmed');
      
      let summary = "";
      
      if (ensemblResult) {
        const additionalInfo = ensemblResult.additionalData;
        summary += `${ensemblResult.title} - ${additionalInfo?.species || 'Species unknown'}${additionalInfo?.chromosome ? `, Chromosome ${additionalInfo.chromosome}` : ''}. `;
      }
      
      if (uniprotResult) {
        summary += `Protein: ${uniprotResult.title} (${uniprotResult.additionalData?.length || 'unknown length'}). `;
      }
      
      if (pdbResult) {
        summary += `Structure available in PDB (Resolution: ${pdbResult.additionalData?.resolution || 'unknown'}). `;
      }
      
      if (pubmedResults.length > 0) {
        summary += `Found ${pubmedResults.length} related publications. `;
      }
      
      if (summary === "" && ncbiResult) {
        summary = `${ncbiResult.title} - ${ncbiResult.additionalData?.organism || ''}. ${ncbiResult.description.split('.')[0]}.`;
      }
      
      return summary || `Found ${results.length} results across multiple biological databases.`;
    }
    
    // Call our Python backend API
    const response = await fetch('http://localhost:5000/api/summary', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ results }),
    });
    
    if (!response.ok) {
      throw new Error(`API error: ${response.status}`);
    }
    
    const data = await response.json();
    return data.summary;
  } catch (error) {
    console.error("Summary error:", error);
    return "Unable to generate summary. Please try again.";
  }
}

// Group results by database source
export function groupResultsBySource(results: ResultItem[]): Record<DatabaseSource, ResultItem[]> {
  const grouped: Record<DatabaseSource, ResultItem[]> = {
    ensembl: [],
    ncbi: [],
    uniprot: [],
    ebi: [],
    ddbj: [],
    github: [],
    pdb: [],
    pubmed: [],
    ucsc: []
  };
  
  results.forEach(result => {
    if (result.source in grouped) {
      grouped[result.source as DatabaseSource].push(result);
    } else {
      console.warn(`Unknown database source: ${result.source}`);
    }
  });
  
  return grouped;
}
