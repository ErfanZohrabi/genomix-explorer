
import React, { useState } from 'react';
import { Header } from '../components/Header';
import { SearchBar } from '../components/SearchBar';
import { EmptyState } from '../components/EmptyState';
import { LoadingState } from '../components/LoadingState';
import { ResultItem, ResultsGrid } from '../components/ResultCard';
import { searchBioDatabases, generateSummary, groupResultsBySource } from '../lib/api';
import { DatabaseSource } from '../components/ResultCard';
import { ChevronDown, ChevronUp, Database, AlertCircle, FileText } from 'lucide-react';
import { Button } from '@/components/ui/button';
import {
  Alert,
  AlertDescription,
  AlertTitle,
} from "@/components/ui/alert";
import { useToast } from "@/hooks/use-toast";

const Index = () => {
  const [query, setQuery] = useState<string>('');
  const [isSearching, setIsSearching] = useState<boolean>(false);
  const [results, setResults] = useState<ResultItem[]>([]);
  const [groupedResults, setGroupedResults] = useState<Record<DatabaseSource, ResultItem[]> | null>(null);
  const [summary, setSummary] = useState<string>('');
  const [expandedSections, setExpandedSections] = useState<Record<string, boolean>>({});
  const [error, setError] = useState<string>('');
  
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
    ensembl: 'Ensembl',
    ncbi: 'NCBI',
    uniprot: 'UniProt',
    ebi: 'EBI ENA',
    ddbj: 'DDBJ',
    github: 'GitHub'
  };

  return (
    <div className="min-h-screen flex flex-col bg-background">
      <Header />
      
      <main className="flex-1 container mx-auto pt-24 px-4 sm:px-6 lg:px-8 max-w-7xl">
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
                <Alert>
                  <FileText className="h-4 w-4" />
                  <AlertTitle>Summary</AlertTitle>
                  <AlertDescription>{summary}</AlertDescription>
                </Alert>
              )}
              
              {groupedResults && Object.entries(groupedResults).map(([source, items]) => {
                if (items.length === 0) return null;
                
                const dbSource = source as DatabaseSource;
                const isExpanded = expandedSections[source] || false;
                
                return (
                  <div key={source} className="rounded-lg border bg-card text-card-foreground shadow-sm">
                    <div 
                      className="flex items-center justify-between p-6 cursor-pointer hover:bg-accent/50 transition-colors"
                      onClick={() => toggleSection(source)}
                    >
                      <h2 className="text-xl font-semibold flex items-center gap-2">
                        <Database size={20} className="text-primary" />
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
              })}
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
