"""
Erowid Experience Vault Scraper
Extracts trip reports for NLP analysis

NOTE: Respect Erowid's terms of service:
- Don't scrape aggressively
- Cache data locally
- Don't redistribute raw data without permission
- Use for research/educational purposes only
"""

import os
import re
import json
import time
import hashlib
from typing import Optional, Dict, List, Generator
from dataclasses import dataclass, asdict
from datetime import datetime
import requests
from bs4 import BeautifulSoup


@dataclass
class ExperienceReport:
    """Structured experience report data"""
    id: str
    title: str
    author: Optional[str]
    substance: str
    substances_detail: List[Dict]  # [{name, amount, route}]
    date_submitted: Optional[str]
    date_experience: Optional[str]
    body_weight: Optional[str]
    gender: Optional[str]
    age: Optional[int]
    year: Optional[int]
    text: str
    categories: List[str]  # e.g., "First Times", "Bad Trips", "Mystical"
    url: str


class ErowidScraper:
    """
    Scraper for Erowid Experience Vaults
    
    Main entry points:
    - /experiences/exp_list.shtml - Browse by substance
    - /experiences/exp.cgi?New - New reports
    - /experiences/subs/exp_*.shtml - Substance-specific pages
    """
    
    BASE_URL = "https://www.erowid.org"
    EXPERIENCES_URL = f"{BASE_URL}/experiences"
    
    # Rate limiting
    REQUEST_DELAY = 2.0  # seconds between requests
    
    def __init__(self, cache_dir: str = "data/raw/erowid/cache"):
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "ConsciousnessResearch/1.0 (Educational/Research)",
            "Accept": "text/html"
        })
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.last_request = 0
    
    def _rate_limit(self):
        """Enforce rate limiting"""
        elapsed = time.time() - self.last_request
        if elapsed < self.REQUEST_DELAY:
            time.sleep(self.REQUEST_DELAY - elapsed)
        self.last_request = time.time()
    
    def _get_cache_path(self, url: str) -> str:
        """Get cache file path for URL"""
        url_hash = hashlib.md5(url.encode()).hexdigest()
        return os.path.join(self.cache_dir, f"{url_hash}.html")
    
    def _fetch(self, url: str, use_cache: bool = True) -> str:
        """Fetch URL with caching"""
        cache_path = self._get_cache_path(url)
        
        # Check cache
        if use_cache and os.path.exists(cache_path):
            with open(cache_path, "r", encoding="utf-8", errors="ignore") as f:
                return f.read()
        
        # Fetch fresh
        self._rate_limit()
        response = self.session.get(url)
        response.raise_for_status()
        
        html = response.text
        
        # Cache result
        with open(cache_path, "w", encoding="utf-8") as f:
            f.write(html)
        
        return html
    
    def get_substance_list(self) -> List[Dict]:
        """Get list of all substances with experience reports"""
        url = f"{self.EXPERIENCES_URL}/exp_substance_list.php"
        html = self._fetch(url)
        soup = BeautifulSoup(html, "html.parser")
        
        substances = []
        # Find substance links in the list
        for link in soup.find_all("a", href=re.compile(r"subs/exp_.*\.shtml")):
            name = link.get_text(strip=True)
            href = link.get("href")
            if name and href:
                substances.append({
                    "name": name,
                    "url": f"{self.EXPERIENCES_URL}/{href}" if not href.startswith("http") else href
                })
        
        return substances
    
    def get_experience_list_for_substance(
        self, 
        substance_url: str,
        max_pages: int = 10
    ) -> List[Dict]:
        """Get list of experience report links for a substance"""
        experiences = []
        page = 1
        
        while page <= max_pages:
            url = f"{substance_url}?Start={(page-1)*25}&Max=25"
            html = self._fetch(url)
            soup = BeautifulSoup(html, "html.parser")
            
            # Find experience links
            found = 0
            for link in soup.find_all("a", href=re.compile(r"exp\.php\?ID=\d+")):
                title = link.get_text(strip=True)
                href = link.get("href")
                
                if title and href:
                    exp_id = re.search(r"ID=(\d+)", href)
                    if exp_id:
                        # Construct proper URL
                        if href.startswith("http"):
                            full_url = href
                        elif href.startswith("/"):
                            full_url = f"{self.BASE_URL}{href}"
                        else:
                            full_url = f"{self.EXPERIENCES_URL}/{href}"
                        experiences.append({
                            "id": exp_id.group(1),
                            "title": title,
                            "url": full_url
                        })
                        found += 1
            
            if found == 0:
                break
            
            page += 1
        
        return experiences
    
    def parse_experience_report(self, url: str) -> Optional[ExperienceReport]:
        """Parse a single experience report page"""
        try:
            html = self._fetch(url)
            soup = BeautifulSoup(html, "html.parser")
            
            # Extract ID from URL
            exp_id = re.search(r"ID=(\d+)", url)
            exp_id = exp_id.group(1) if exp_id else hashlib.md5(url.encode()).hexdigest()[:8]
            
            # Title
            title_elem = soup.find("div", class_="title")
            title = title_elem.get_text(strip=True) if title_elem else "Unknown"
            
            # Author
            author_elem = soup.find("div", class_="author")
            author = author_elem.get_text(strip=True).replace("by ", "") if author_elem else None
            
            # Substance (from heading or breadcrumbs)
            substance = "Unknown"
            substance_elem = soup.find("div", class_="substance")
            if substance_elem:
                substance = substance_elem.get_text(strip=True)
            
            # Detailed substance info from the dose table
            substances_detail = []
            dose_table = soup.find("table", class_="dosechart")
            if dose_table:
                for row in dose_table.find_all("tr"):
                    cols = row.find_all("td")
                    if len(cols) >= 3:
                        substances_detail.append({
                            "time": cols[0].get_text(strip=True) if len(cols) > 0 else None,
                            "amount": cols[1].get_text(strip=True) if len(cols) > 1 else None,
                            "substance": cols[2].get_text(strip=True) if len(cols) > 2 else None,
                            "route": cols[3].get_text(strip=True) if len(cols) > 3 else None
                        })
            
            # Body weight
            body_weight = None
            weight_elem = soup.find(string=re.compile(r"Body weight"))
            if weight_elem:
                weight_match = re.search(r"(\d+)\s*(kg|lb)", weight_elem)
                if weight_match:
                    body_weight = f"{weight_match.group(1)} {weight_match.group(2)}"
            
            # Main text
            text = ""
            report_text = soup.find("div", class_="report-text-surround")
            if report_text:
                text = report_text.get_text(separator="\n", strip=True)
            
            # Categories/tags
            categories = []
            for cat_link in soup.find_all("a", href=re.compile(r"exp_list\.shtml\?.*Category")):
                cat_name = cat_link.get_text(strip=True)
                if cat_name:
                    categories.append(cat_name)
            
            return ExperienceReport(
                id=exp_id,
                title=title,
                author=author,
                substance=substance,
                substances_detail=substances_detail,
                date_submitted=None,  # Would need additional parsing
                date_experience=None,
                body_weight=body_weight,
                gender=None,
                age=None,
                year=None,
                text=text,
                categories=categories,
                url=url
            )
            
        except Exception as e:
            print(f"Error parsing {url}: {e}")
            return None
    
    def scrape_substance(
        self,
        substance_name: str,
        substance_url: str,
        max_reports: int = 100,
        output_dir: str = "data/raw/erowid"
    ) -> List[ExperienceReport]:
        """Scrape all experience reports for a substance"""
        os.makedirs(f"{output_dir}/reports", exist_ok=True)
        
        print(f"Fetching report list for {substance_name}...")
        report_list = self.get_experience_list_for_substance(substance_url)
        print(f"Found {len(report_list)} reports")
        
        reports = []
        for i, item in enumerate(report_list[:max_reports]):
            print(f"[{i+1}/{min(len(report_list), max_reports)}] {item['title'][:50]}...")
            
            report = self.parse_experience_report(item["url"])
            if report:
                reports.append(report)
                
                # Save individual report
                with open(f"{output_dir}/reports/{report.id}.json", "w") as f:
                    json.dump(asdict(report), f, indent=2)
        
        # Save summary
        safe_name = substance_name.replace("/", "_").replace(" ", "_")
        with open(f"{output_dir}/{safe_name}_reports.json", "w") as f:
            json.dump([asdict(r) for r in reports], f, indent=2)
        
        return reports


def create_corpus_index(reports_dir: str, output_file: str):
    """Create searchable index of all reports"""
    index = []
    
    for filename in os.listdir(reports_dir):
        if filename.endswith(".json"):
            with open(os.path.join(reports_dir, filename)) as f:
                report = json.load(f)
                index.append({
                    "id": report["id"],
                    "title": report["title"],
                    "substance": report["substance"],
                    "substances_detail": report.get("substances_detail", []),
                    "categories": report.get("categories", []),
                    "text_length": len(report.get("text", "")),
                    "file": filename
                })
    
    with open(output_file, "w") as f:
        json.dump(index, f, indent=2)
    
    print(f"Indexed {len(index)} reports to {output_file}")
    return index


if __name__ == "__main__":
    scraper = ErowidScraper()
    
    # Example: Get substance list
    print("Fetching substance list...")
    substances = scraper.get_substance_list()
    print(f"Found {len(substances)} substances")
    for s in substances[:10]:
        print(f"  - {s['name']}")
    
    # Example: Parse a single report (using a generic URL pattern)
    # Note: Replace with actual report URL for testing
    # report = scraper.parse_experience_report("https://www.erowid.org/experiences/exp.php?ID=12345")
    # if report:
    #     print(json.dumps(asdict(report), indent=2))
    
    # CAUTION: Only run full scraping with permission and appropriate delays
    # scraper.scrape_substance("LSD", "https://www.erowid.org/experiences/subs/exp_LSD.shtml")
