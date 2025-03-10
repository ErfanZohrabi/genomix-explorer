
import React from 'react';
import { Link } from 'react-router-dom';
import { ThemeToggle } from './ThemeToggle';
import { BookOpen, HelpCircle, Home } from 'lucide-react';
import { Button } from "@/components/ui/button";

export function Header() {
  return (
    <header className="w-full py-4 px-6 md:px-8 flex items-center justify-between border-b border-border/50 bg-background/80 backdrop-blur-md fixed top-0 left-0 right-0 z-50 transition-all duration-300">
      <div className="flex items-center gap-2">
        <Link to="/" className="flex items-center gap-2">
          <span className="w-9 h-9 rounded-lg bg-biosearch-600 flex items-center justify-center">
            <span className="text-white font-semibold text-lg">B</span>
          </span>
          <span className="font-semibold text-lg hidden md:block">Bio_Search</span>
        </Link>
      </div>
      
      <div className="flex items-center gap-2">
        <Link to="/">
          <Button variant="ghost" size="sm" className="hidden md:flex items-center gap-1">
            <Home size={16} />
            <span>Home</span>
          </Button>
        </Link>
        <Link to="/about">
          <Button variant="ghost" size="sm" className="hidden md:flex items-center gap-1">
            <BookOpen size={16} />
            <span>About</span>
          </Button>
        </Link>
        <Link to="/">
          <Button variant="ghost" size="icon" className="md:hidden rounded-full w-9 h-9">
            <Home size={18} />
          </Button>
        </Link>
        <Link to="/about">
          <Button variant="ghost" size="icon" className="md:hidden rounded-full w-9 h-9">
            <HelpCircle size={18} />
          </Button>
        </Link>
        <ThemeToggle />
      </div>
    </header>
  );
}
