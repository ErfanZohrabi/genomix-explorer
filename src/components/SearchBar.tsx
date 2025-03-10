
import React, { useState, useRef, useEffect } from 'react';
import { Search, X } from 'lucide-react';
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";

interface SearchBarProps {
  onSearch: (query: string) => void;
  isSearching: boolean;
  className?: string;
}

export function SearchBar({ onSearch, isSearching, className = '' }: SearchBarProps) {
  const [query, setQuery] = useState('');
  const [hasFocus, setHasFocus] = useState(false);
  const inputRef = useRef<HTMLInputElement>(null);
  const searchTimeout = useRef<NodeJS.Timeout | null>(null);

  const handleSearch = () => {
    if (query.trim()) {
      onSearch(query.trim());
    }
  };

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      handleSearch();
    }
  };

  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newQuery = e.target.value;
    setQuery(newQuery);
    
    // Debounce search for suggestions
    if (searchTimeout.current) {
      clearTimeout(searchTimeout.current);
    }
    
    // Only trigger for non-empty queries with at least 3 characters
    if (newQuery.trim().length >= 3) {
      searchTimeout.current = setTimeout(() => {
        // This could trigger a suggestion search in the future
      }, 300);
    }
  };

  const clearSearch = () => {
    setQuery('');
    if (inputRef.current) {
      inputRef.current.focus();
    }
  };

  useEffect(() => {
    // Clean up timeout on unmount
    return () => {
      if (searchTimeout.current) {
        clearTimeout(searchTimeout.current);
      }
    };
  }, []);

  return (
    <div className={`relative ${className}`}>
      <div 
        className={`flex items-center bg-card border transition-all duration-300 rounded-full overflow-hidden shadow-sm 
        ${hasFocus ? 'ring-2 ring-primary/20 border-primary/30' : 'border-border hover:border-input'}`}
      >
        <div className="pl-4">
          <Search size={20} className="text-muted-foreground" />
        </div>
        <Input
          ref={inputRef}
          type="text"
          placeholder="Search genes, proteins, or sequences..."
          value={query}
          onChange={handleChange}
          onKeyDown={handleKeyDown}
          onFocus={() => setHasFocus(true)}
          onBlur={() => setHasFocus(false)}
          className="border-0 shadow-none focus-visible:ring-0 bg-transparent rounded-none h-12 px-3 text-base"
        />
        {query && (
          <Button 
            variant="ghost" 
            size="icon" 
            onClick={clearSearch} 
            className="rounded-full mr-1 hover:bg-accent/50"
          >
            <X size={16} className="text-muted-foreground" />
            <span className="sr-only">Clear search</span>
          </Button>
        )}
        <Button 
          onClick={handleSearch}
          disabled={!query.trim() || isSearching}
          variant="default"
          className="m-1 h-10 px-4 rounded-full font-medium"
        >
          {isSearching ? 'Searching...' : 'Search'}
        </Button>
      </div>
    </div>
  );
}
