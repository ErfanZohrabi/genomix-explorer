
import React, { useState } from 'react';
import { Header } from '../components/Header';
import { SearchBar } from '../components/SearchBar';
import { EmptyState } from '../components/EmptyState';
import { LoadingState } from '../components/LoadingState';
import { ResultItem, ResultsGrid } from '../components/ResultCard';
import { searchBioDatabases, generateSummary, groupResultsBySource } from '../lib/api';
import { DatabaseSource } from '../components/ResultCard';
import { ChevronDown, ChevronUp, Database, AlertCircle, FileText, Search, Dna, Github, FileCode, Book, BookOpen } from 'lucide-react';
import { Button } from '@/components/ui/button';
import {
  Alert,
  AlertDescription,
  AlertTitle,
} from "@/components/ui/alert";
import { useToast } from "@/hooks/use-toast";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";

const Index = () => {
  const [query, setQuery] = useState<string>('');
  const [isSearching, setIsSearching] = useState<boolean>(false);
  const [results, setResults] = useState<ResultItem[]>([]);
  const [groupedResults, setGroupedResults] = useState<Record<DatabaseSource, ResultItem[]> | null>(null);
  const [summary, setSummary] = useState<string>('');
  const [expandedSections, setExpandedSections] = useState<Record<string, boolean>>({});
  const [error, setError] = useState<string>('');
  const [activeTab, setActiveTab] = useState<string>("all");
  
  const { toast } = useToast();
  
  const handleSearch = async (searchQuery: string) => {
    setQuery(searchQuery);
    setIsSearching(true);
    setError('');
    
    try {
      // Show a toast notification that the search has started
      toast({
        title: "Searching databases...",
        description: "This might take a few seconds.",
      });
      
      const searchResults = await searchBioDatabases(searchQuery);
      setResults(searchResults);
      
      const grouped = groupResultsBySource(searchResults);
      setGroupedResults(grouped);
      
      // Generate and set the summary
      const summaryText = await generateSummary(searchResults);
      setSummary(summaryText);
      
      // Initialize all non-empty sections as expanded
      const initialExpanded: Record<string, boolean> = {};
      Object.entries(grouped).forEach(([source, items]) => {
        initialExpanded[source] = items.length > 0;
      });
      setExpandedSections(initialExpanded);
      
      // Show success toast if results were found
      if (searchResults.length > 0) {
        toast({
          title: "Search complete",
          description: `Found ${searchResults.length} results across multiple databases.`,
        });
      } else {
        toast({
          variant: "destructive",
          title: "No results found",
          description: "Try refining your search terms.",
        });
      }
    } catch (error) {
      console.error('Search error:', error);
      setResults([]);
      setGroupedResults(null);
      setSummary('');
      setError('An error occurred while searching. Please try again.');
      
      toast({
        variant: "destructive",
        title: "Search failed",
        description: "An error occurred while searching the databases.",
      });
    } finally {
      setIsSearching(false);
    }
  };
  
  const toggleSection = (section: string) => {
    setExpandedSections(prev => ({
      ...prev,
      [section]: !prev[section]
    }));
  };
  
  const databaseLabels: Record<DatabaseSource, string> = {
    ensembl: 'Ensembl Gene Database',
    ncbi: 'NCBI Gene & Nucleotide',
    uniprot: 'UniProt Protein Database',
    ebi: 'EBI ENA Sequence Database',
    ddbj: 'DDBJ Nucleotide Database',
    github: 'GitHub Repositories',
    pdb: 'PDB Protein Structures',
    pubmed: 'PubMed Literature',
    ucsc: 'UCSC Genome Browser'
  };
  
  const databaseIcons: Record<DatabaseSource, React.ReactNode> = {
    ensembl: <Dna size={20} className="text-blue-500" />,
    ncbi: <Database size={20} className="text-green-600" />,
    uniprot: <FileCode size={20} className="text-yellow-600" />,
    ebi: <FileText size={20} className="text-purple-600" />,
    ddbj: <Dna size={20} className="text-indigo-600" />,
    github: <Github size={20} className="text-gray-600" />,
    pdb: <FileCode size={20} className="text-pink-600" />,
    pubmed: <BookOpen size={20} className="text-red-600" />,
    ucsc: <Search size={20} className="text-teal-600" />
  };
  
  const getResultCount = () => {
    if (!groupedResults) return 0;
    return Object.values(groupedResults).reduce((total, items) => total + items.length, 0);
  };
  
  const getTotalResultsPerTab = (tab: string): number => {
    if (!groupedResults) return 0;
    if (tab === "all") return getResultCount();
    
    if (tab === "genes") {
      return (groupedResults.ensembl?.length || 0) + (groupedResults.ncbi?.length || 0);
    }
    if (tab === "proteins") {
      return (groupedResults.uniprot?.length || 0) + (groupedResults.pdb?.length || 0);
    }
    if (tab === "sequences") {
      return (groupedResults.ebi?.length || 0) + (groupedResults.ddbj?.length || 0);
    }
    if (tab === "papers") {
      return groupedResults.pubmed?.length || 0;
    }
    if (tab === "repos") {
      return groupedResults.github?.length || 0;
    }
    
    return 0;
  };
  
  const renderTabContent = (tab: string) => {
    if (!groupedResults) return null;
    
    const filteredSources = Object.entries(groupedResults).filter(([source, items]) => {
      if (items.length === 0) return false;
      
      if (tab === "all") return true;
      if (tab === "genes" && (source === "ensembl" || source === "ncbi")) return true;
      if (tab === "proteins" && (source === "uniprot" || source === "pdb")) return true;
      if (tab === "sequences" && (source === "ebi" || source === "ddbj")) return true;
      if (tab === "papers" && source === "pubmed") return true;
      if (tab === "repos" && source === "github") return true;
      
      return false;
    });
    
    return filteredSources.map(([source, items]) => {
      const dbSource = source as DatabaseSource;
      const isExpanded = expandedSections[source] || false;
      
      return (
        <div key={source} className="rounded-lg border bg-card text-card-foreground shadow-sm mb-6">
          <div 
            className="flex items-center justify-between p-6 cursor-pointer hover:bg-accent/50 transition-colors"
            onClick={() => toggleSection(source)}
          >
            <h2 className="text-xl font-semibold flex items-center gap-2">
              {databaseIcons[dbSource]}
              <span>{databaseLabels[dbSource]} Results</span>
              <span className="ml-2 text-sm font-normal text-muted-foreground">
                ({items.length} {items.length === 1 ? 'item' : 'items'})
              </span>
            </h2>
            <Button variant="ghost" size="sm" className="h-8 w-8 p-0">
              {isExpanded ? <ChevronUp size={18} /> : <ChevronDown size={18} />}
            </Button>
          </div>
          
          {isExpanded && (
            <div className="p-6 pt-0">
              <ResultsGrid items={items} />
            </div>
          )}
        </div>
      );
    });
  };

  return (
    <div className="min-h-screen flex flex-col bg-background">
      <Header />
      
      <main className="flex-1 container mx-auto pt-16 px-4 sm:px-6 lg:px-8 max-w-7xl">
        <div className="text-center mb-8">
          <h1 className="text-4xl font-bold mb-3">BioSearch</h1>
          <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
            Search across multiple biological databases, including genes, proteins, sequences, structures, and literature.
          </p>
        </div>
        
        <div className={`search-container transition-all duration-500 ease-out ${results.length > 0 ? 'searching' : ''}`}>
          <SearchBar 
            onSearch={handleSearch} 
            isSearching={isSearching} 
            className="mb-8"
          />
          
          {error && (
            <Alert variant="destructive" className="mb-8">
              <AlertCircle className="h-4 w-4" />
              <AlertTitle>Error</AlertTitle>
              <AlertDescription>{error}</AlertDescription>
            </Alert>
          )}
          
          {isSearching ? (
            <LoadingState />
          ) : results.length > 0 ? (
            <div className="animate-fade-in space-y-8 pb-16">
              {summary && (
                <Alert className="mb-8 bg-blue-50 dark:bg-blue-900/30 border-blue-300 dark:border-blue-800">
                  <FileText className="h-4 w-4 text-blue-600 dark:text-blue-400" />
                  <AlertTitle className="text-blue-800 dark:text-blue-300">Summary</AlertTitle>
                  <AlertDescription className="text-blue-700 dark:text-blue-400">{summary}</AlertDescription>
                </Alert>
              )}
              
              <Tabs defaultValue="all" value={activeTab} onValueChange={setActiveTab} className="w-full">
                <div className="flex justify-center mb-8">
                  <TabsList className="grid grid-cols-2 md:grid-cols-6 max-w-xl">
                    <TabsTrigger value="all" className="flex gap-1 items-center">
                      <Database size={16} />
                      <span>All ({getTotalResultsPerTab('all')})</span>
                    </TabsTrigger>
                    <TabsTrigger value="genes" className="flex gap-1 items-center">
                      <Dna size={16} />
                      <span>Genes ({getTotalResultsPerTab('genes')})</span>
                    </TabsTrigger>
                    <TabsTrigger value="proteins" className="flex gap-1 items-center">
                      <FileCode size={16} />
                      <span>Proteins ({getTotalResultsPerTab('proteins')})</span>
                    </TabsTrigger>
                    <TabsTrigger value="sequences" className="flex gap-1 items-center">
                      <FileText size={16} />
                      <span>Sequences ({getTotalResultsPerTab('sequences')})</span>
                    </TabsTrigger>
                    <TabsTrigger value="papers" className="flex gap-1 items-center">
                      <Book size={16} />
                      <span>Papers ({getTotalResultsPerTab('papers')})</span>
                    </TabsTrigger>
                    <TabsTrigger value="repos" className="flex gap-1 items-center">
                      <Github size={16} />
                      <span>GitHub ({getTotalResultsPerTab('repos')})</span>
                    </TabsTrigger>
                  </TabsList>
                </div>
                
                <TabsContent value="all">
                  {renderTabContent("all")}
                </TabsContent>
                
                <TabsContent value="genes">
                  {renderTabContent("genes")}
                </TabsContent>
                
                <TabsContent value="proteins">
                  {renderTabContent("proteins")}
                </TabsContent>
                
                <TabsContent value="sequences">
                  {renderTabContent("sequences")}
                </TabsContent>
                
                <TabsContent value="papers">
                  {renderTabContent("papers")}
                </TabsContent>
                
                <TabsContent value="repos">
                  {renderTabContent("repos")}
                </TabsContent>
              </Tabs>
            </div>
          ) : (
            <EmptyState />
          )}
        </div>
      </main>
    </div>
  );
};

export default Index;
