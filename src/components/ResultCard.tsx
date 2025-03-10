
import React from 'react';
import { ExternalLink } from 'lucide-react';
import { Card, CardContent } from "@/components/ui/card";

// Define types for different database sources
export type DatabaseSource = 'ensembl' | 'ncbi' | 'uniprot' | 'ebi' | 'ddbj' | 'github';

export interface ResultItem {
  id: string;
  title: string;
  description: string;
  source: DatabaseSource;
  url: string;
  additionalData?: Record<string, any>;
}

interface ResultCardProps {
  item: ResultItem;
  index: number;
}

// Map of icons/colors for different sources
const sourceConfig: Record<DatabaseSource, { color: string; label: string }> = {
  ensembl: { color: 'bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-300', label: 'Ensembl' },
  ncbi: { color: 'bg-amber-100 text-amber-800 dark:bg-amber-900 dark:text-amber-300', label: 'NCBI' },
  uniprot: { color: 'bg-purple-100 text-purple-800 dark:bg-purple-900 dark:text-purple-300', label: 'UniProt' },
  ebi: { color: 'bg-emerald-100 text-emerald-800 dark:bg-emerald-900 dark:text-emerald-300', label: 'EBI ENA' },
  ddbj: { color: 'bg-rose-100 text-rose-800 dark:bg-rose-900 dark:text-rose-300', label: 'DDBJ' },
  github: { color: 'bg-gray-100 text-gray-800 dark:bg-gray-800 dark:text-gray-300', label: 'GitHub' }
};

export function ResultCard({ item, index }: ResultCardProps) {
  const sourceStyle = sourceConfig[item.source];
  
  return (
    <div 
      className="animate-slide-up" 
      style={{ animationDelay: `${index * 0.05}s` }}
    >
      <Card className="glass-card overflow-hidden h-full">
        <CardContent className="p-4 flex flex-col h-full">
          <div className="flex justify-between items-start gap-2 mb-2">
            <span className={`px-2 py-1 text-xs font-medium rounded-md ${sourceStyle.color}`}>
              {sourceStyle.label}
            </span>
            <a 
              href={item.url} 
              target="_blank" 
              rel="noopener noreferrer"
              className="text-muted-foreground hover:text-foreground transition-colors"
              aria-label={`Open ${item.title} in ${sourceStyle.label}`}
            >
              <ExternalLink size={16} />
            </a>
          </div>
          
          <h3 className="font-medium text-lg mb-1 line-clamp-1">{item.title}</h3>
          <p className="text-muted-foreground text-sm flex-grow line-clamp-3 mb-2">{item.description}</p>
          
          {item.additionalData && Object.keys(item.additionalData).length > 0 && (
            <div className="mt-2 pt-2 border-t border-border/50 text-xs text-muted-foreground">
              {Object.entries(item.additionalData).map(([key, value]) => (
                <div key={key} className="flex items-center justify-between">
                  <span className="font-medium capitalize">{key}:</span> 
                  <span>{String(value)}</span>
                </div>
              ))}
            </div>
          )}
        </CardContent>
      </Card>
    </div>
  );
}

export function ResultsGrid({ items }: { items: ResultItem[] }) {
  return (
    <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 w-full">
      {items.map((item, index) => (
        <ResultCard key={item.id} item={item} index={index} />
      ))}
    </div>
  );
}
