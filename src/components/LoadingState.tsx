
import React from 'react';

export function LoadingState() {
  return (
    <div className="w-full py-16 animate-fade-in">
      <div className="flex flex-col items-center justify-center gap-4">
        <div className="relative w-16 h-16">
          <div className="absolute inset-0 border-t-4 border-biosearch-500 border-solid rounded-full animate-spin"></div>
          <div className="absolute inset-2 border-r-4 border-biosearch-300 border-solid rounded-full animate-spin animate-[spin_1.5s_linear_infinite]"></div>
        </div>
        <div className="text-lg font-medium text-center">Searching databases...</div>
        <div className="text-sm text-muted-foreground text-center max-w-sm">
          Querying Ensembl, NCBI, UniProt, EBI ENA, DDBJ, and GitHub repositories
        </div>
        <div className="text-xs text-muted-foreground text-center max-w-xs mt-2">
          Each database API is being queried in parallel to bring you comprehensive results
        </div>
      </div>
    </div>
  );
}
