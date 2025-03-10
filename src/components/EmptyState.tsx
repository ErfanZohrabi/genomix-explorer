
import React from 'react';
import { Database, Dna, FileSearch, Github } from 'lucide-react';

export function EmptyState() {
  return (
    <div className="flex flex-col items-center justify-center pt-10 pb-16 animate-fade-in">
      <div className="text-center max-w-2xl space-y-6">
        <h1 className="text-2xl sm:text-3xl font-semibold tracking-tight mb-3">
          Discover Biological Data Across Multiple Databases
        </h1>
        <p className="text-muted-foreground text-base sm:text-lg">
          Search across Ensembl, NCBI, UniProt, EBI ENA, DDBJ, and GitHub repositories
          to find genes, proteins, and biological sequences.
        </p>
        <div className="grid grid-cols-2 md:grid-cols-3 gap-4 pt-6">
          <Database className="grid-database mx-auto text-biosearch-500 w-12 h-12 opacity-90" />
          <Dna className="grid-dna mx-auto text-biosearch-600 w-12 h-12 opacity-90" />
          <FileSearch className="grid-search mx-auto text-biosearch-700 w-12 h-12 opacity-90" />
          <div className="hidden md:block"></div>
          <Github className="grid-github mx-auto text-biosearch-800 w-12 h-12 opacity-90" />
          <div className="hidden md:block"></div>
        </div>
        <div className="pt-6">
          <h2 className="text-xl font-medium mb-3">Try searching for:</h2>
          <div className="flex flex-wrap justify-center gap-2">
            <ExampleTag text="BRCA1" />
            <ExampleTag text="p53" />
            <ExampleTag text="insulin" />
            <ExampleTag text="ATGGCGCGAT..." />
            <ExampleTag text="coronavirus" />
          </div>
        </div>
      </div>
    </div>
  );
}

function ExampleTag({ text }: { text: string }) {
  return (
    <div className="px-3 py-1.5 bg-secondary/70 hover:bg-secondary text-secondary-foreground rounded-full text-sm transition-colors cursor-pointer">
      {text}
    </div>
  );
}
