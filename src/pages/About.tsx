
import React from 'react';
import { Header } from '../components/Header';
import { Database, Github, Globe, BookOpen, Code, Server } from 'lucide-react';

const About = () => {
  return (
    <div className="min-h-screen flex flex-col bg-background">
      <Header />
      
      <main className="flex-1 pt-24 pb-16 px-4 sm:px-6 lg:px-8 max-w-6xl mx-auto">
        <div className="mx-auto max-w-4xl animate-fade-in">
          <div className="space-y-8">
            <div className="text-center mb-10">
              <h1 className="text-3xl md:text-4xl font-bold tracking-tight mb-4">About Bio_Search</h1>
              <p className="text-lg text-muted-foreground">
                A powerful platform to explore bioinformatics data across multiple databases
              </p>
            </div>
            
            <section className="glass-card p-6 md:p-8 rounded-2xl mb-8">
              <h2 className="text-2xl font-semibold mb-4 flex items-center gap-2">
                <Globe className="text-biosearch-600" size={24} />
                <span>Our Mission</span>
              </h2>
              <p className="text-muted-foreground mb-4">
                Bio_Search aims to simplify the process of finding and analyzing biological data by providing
                a unified search interface across multiple authoritative bioinformatics databases.
              </p>
              <p className="text-muted-foreground">
                Whether you're researching genes, proteins, or biological sequences, our platform helps you
                discover relevant information quickly and efficiently from sources like Ensembl, NCBI, UniProt,
                EBI ENA, DDBJ, and GitHub repositories.
              </p>
            </section>
            
            <section className="glass-card p-6 md:p-8 rounded-2xl mb-8">
              <h2 className="text-2xl font-semibold mb-6 flex items-center gap-2">
                <Database className="text-biosearch-600" size={24} />
                <span>Supported Databases</span>
              </h2>
              
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                <DatabaseCard 
                  name="Ensembl" 
                  description="Genome browser and annotation system focusing on vertebrate genomes"
                  url="https://ensembl.org/"
                />
                <DatabaseCard 
                  name="NCBI" 
                  description="National Center for Biotechnology Information providing access to biomedical and genomic information"
                  url="https://www.ncbi.nlm.nih.gov/"
                />
                <DatabaseCard 
                  name="UniProt" 
                  description="Comprehensive resource for protein sequence and functional information"
                  url="https://www.uniprot.org/"
                />
                <DatabaseCard 
                  name="EBI ENA" 
                  description="European Nucleotide Archive providing access to nucleotide sequencing information"
                  url="https://www.ebi.ac.uk/ena/browser/home"
                />
                <DatabaseCard 
                  name="DDBJ" 
                  description="DNA Data Bank of Japan collecting nucleotide sequence data"
                  url="https://www.ddbj.nig.ac.jp/"
                />
                <DatabaseCard 
                  name="GitHub" 
                  description="Repositories containing bioinformatics tools, pipelines, and analyses"
                  url="https://github.com/topics/bioinformatics"
                />
              </div>
            </section>
            
            <section className="glass-card p-6 md:p-8 rounded-2xl mb-8">
              <h2 className="text-2xl font-semibold mb-4 flex items-center gap-2">
                <Code className="text-biosearch-600" size={24} />
                <span>Technology Stack</span>
              </h2>
              <p className="text-muted-foreground mb-6">
                Bio_Search is built using modern web technologies to ensure a fast, responsive, and 
                user-friendly experience:
              </p>
              
              <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
                <TechCard name="React" description="Component-based UI library" />
                <TechCard name="TypeScript" description="Type-safe JavaScript" />
                <TechCard name="TailwindCSS" description="Utility-first CSS framework" />
                <TechCard name="React Query" description="Data synchronization library" />
                <TechCard name="ShadCN UI" description="Accessible UI components" />
                <TechCard name="Vite" description="Fast frontend build tool" />
              </div>
            </section>
            
            <section className="glass-card p-6 md:p-8 rounded-2xl">
              <h2 className="text-2xl font-semibold mb-4 flex items-center gap-2">
                <BookOpen className="text-biosearch-600" size={24} />
                <span>How to Use</span>
              </h2>
              <p className="text-muted-foreground mb-6">
                Using Bio_Search is straightforward:
              </p>
              
              <ol className="space-y-4 text-muted-foreground list-decimal pl-5">
                <li>
                  <span className="font-medium text-foreground">Enter your query</span> in the search bar
                  (e.g., gene name, protein identifier, or sequence).
                </li>
                <li>
                  <span className="font-medium text-foreground">Review the results</span> from multiple databases,
                  organized by source.
                </li>
                <li>
                  <span className="font-medium text-foreground">Click on result cards</span> to access the original
                  entries in their respective databases.
                </li>
                <li>
                  <span className="font-medium text-foreground">Read the summary</span> for a quick overview of
                  the most relevant information.
                </li>
              </ol>
            </section>
          </div>
        </div>
      </main>
      
      <footer className="border-t border-border/50 py-6 px-4 text-center text-sm text-muted-foreground">
        <div className="max-w-6xl mx-auto flex flex-col md:flex-row justify-between items-center gap-4">
          <div>© 2023 Bio_Search. All rights reserved.</div>
          <div className="flex items-center gap-4">
            <a href="#" className="hover:text-foreground transition-colors">Terms</a>
            <a href="#" className="hover:text-foreground transition-colors">Privacy</a>
            <a href="https://github.com" target="_blank" rel="noopener noreferrer" className="hover:text-foreground transition-colors">
              <Github size={16} />
            </a>
          </div>
        </div>
      </footer>
    </div>
  );
};

function DatabaseCard({ name, description, url }: { name: string; description: string; url: string }) {
  return (
    <a 
      href={url} 
      target="_blank" 
      rel="noopener noreferrer"
      className="block p-4 rounded-xl border border-border bg-card/50 hover:bg-card transition-colors"
    >
      <h3 className="font-medium text-lg mb-1">{name}</h3>
      <p className="text-sm text-muted-foreground">{description}</p>
      <div className="mt-2 text-xs text-biosearch-600 dark:text-biosearch-400 flex items-center gap-1">
        <span>Visit database</span>
        <Server size={12} />
      </div>
    </a>
  );
}

function TechCard({ name, description }: { name: string; description: string }) {
  return (
    <div className="p-3 rounded-lg border border-border/70 bg-card/30">
      <h3 className="font-medium">{name}</h3>
      <p className="text-xs text-muted-foreground">{description}</p>
    </div>
  );
}

export default About;
