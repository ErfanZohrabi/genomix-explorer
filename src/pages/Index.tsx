
import React, { useState } from 'react';
import { Header } from '../components/Header';
import { SearchBar } from '../components/SearchBar';
import { EmptyState } from '../components/EmptyState';
import { LoadingState } from '../components/LoadingState';
import { ResultItem, ResultsGrid } from '../components/ResultCard';
import { searchBioDatabases, generateSummary, groupResultsBySource } from '../lib/api';
import { DatabaseSource } from '../components/ResultCard';
import { ChevronDown, ChevronUp, Database } from 'lucide-react';
import { Button } from '@/components/ui/button';

const Index = () => {
  const [query, setQuery] = useState<string>('');
  const [isSearching, setIsSearching] = useState<boolean>(false);
  const [results, setResults] = useState<ResultItem[]>([]);
  const [groupedResults, setGroupedResults] = useState<Record<DatabaseSource, ResultItem[]> | null>(null);
  const [summary, setSummary] = useState<string>('');
  const [expandedSections, setExpandedSections] = useState<Record<string, boolean>>({});
  
  const handleSearch = async (searchQuery: string) => {
    setQuery(searchQuery);
    setIsSearching(true);
    
    try {
      const searchResults = await searchBioDatabases(searchQuery);
      setResults(searchResults);
      
      const grouped = groupResultsBySource(searchResults);
      setGroupedResults(grouped);
      
      const summaryText = generateSummary(searchResults);
      setSummary(summaryText);
      
      // Initialize all non-empty sections as expanded
      const initialExpanded: Record<string, boolean> = {};
      Object.entries(grouped).forEach(([source, items]) => {
        initialExpanded[source] = items.length > 0;
      });
      setExpandedSections(initialExpanded);
    } catch (error) {
      console.error('Search error:', error);
      setResults([]);
      setGroupedResults(null);
      setSummary('An error occurred while searching. Please try again.');
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
      
      <main className="flex-1 pt-24 px-4 sm:px-6 lg:px-8">
        <div className={`search-container transition-all duration-500 ease-out ${results.length > 0 ? 'searching' : ''}`}>
          <SearchBar 
            onSearch={handleSearch} 
            isSearching={isSearching} 
            className="mb-8"
          />
          
          {isSearching ? (
            <LoadingState />
          ) : results.length > 0 ? (
            <div className="animate-fade-in space-y-8 pb-16">
              {summary && (
                <div className="glass-card p-4 rounded-xl mb-8">
                  <h2 className="text-lg font-medium mb-1">Summary</h2>
                  <p className="text-muted-foreground">{summary}</p>
                </div>
              )}
              
              {groupedResults && Object.entries(groupedResults).map(([source, items]) => {
                if (items.length === 0) return null;
                
                const dbSource = source as DatabaseSource;
                const isExpanded = expandedSections[source] || false;
                
                return (
                  <div key={source} className="mb-8">
                    <div 
                      className="flex items-center justify-between mb-4 cursor-pointer"
                      onClick={() => toggleSection(source)}
                    >
                      <h2 className="text-xl font-semibold flex items-center gap-2">
                        <Database size={20} className="text-biosearch-600" />
                        <span>{databaseLabels[dbSource]} Results</span>
                        <span className="ml-2 text-sm font-normal text-muted-foreground">
                          ({items.length} {items.length === 1 ? 'item' : 'items'})
                        </span>
                      </h2>
                      <Button variant="ghost" size="sm" className="h-8 w-8 p-0 rounded-full">
                        {isExpanded ? <ChevronUp size={18} /> : <ChevronDown size={18} />}
                      </Button>
                    </div>
                    
                    {isExpanded && (
                      <ResultsGrid items={items} />
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
