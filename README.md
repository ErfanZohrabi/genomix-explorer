
# Bio_Search

A bioinformatics search platform that allows users to search for genes, proteins, and sequences across multiple databases (Ensembl, NCBI, UniProt, EBI ENA, DDBJ, and GitHub).

## Features

- Search across multiple biological databases
- View results organized by database source
- Get concise summaries of search results
- Dark/light mode toggle
- Responsive design

## Tech Stack

- Frontend: React, TypeScript, TailwindCSS
- Backend: Python Flask
- API Integrations: Ensembl, NCBI, UniProt, EBI ENA, DDBJ, GitHub

## Setup Instructions

### Backend (Python Flask)

1. Navigate to the `api` directory:
   ```
   cd api
   ```

2. Create a virtual environment (optional but recommended):
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install dependencies:
   ```
   pip install -r requirements.txt
   ```

4. Run the Flask server:
   ```
   python app.py
   ```
   The server will start on http://localhost:5000

### Frontend (React)

1. In the project root directory, install dependencies:
   ```
   npm install
   ```

2. Run the development server:
   ```
   npm run dev
   ```
   The frontend will start on http://localhost:5173

## Usage

1. Enter a search term in the search bar (e.g., "BRCA1", "p53", "insulin")
2. View the results organized by database source
3. Check the summary for key information about your search term

## API Endpoints

- GET `/api/search?query=<search_term>` - Search across all databases
- POST `/api/summary` - Generate a summary from search results
