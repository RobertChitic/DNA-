#!/usr/bin/env python3
"""
Script: 10_generate_report.py
Description: Generate interactive HTML dashboard using Plotly

Creates:
- Interactive PCA scatter plot
- ADMIXTURE component visualization
- Country probability charts
- Downloadable results summary

Usage:
    python 10_generate_report.py --results-dir <dir> --output <file.html>
"""

import argparse
import os
import sys
import logging
from typing import Dict, List, Optional
from datetime import datetime

import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
except ImportError:
    print("Plotly not installed. Install with: pip install plotly")
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# HTML template
HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Ancestry Analysis Results</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background-color: #f5f5f5;
            color: #333;
            line-height: 1.6;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px 20px;
            text-align: center;
        }}
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        .header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        .section {{
            background: white;
            border-radius: 10px;
            padding: 30px;
            margin: 20px 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .section h2 {{
            color: #667eea;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #eee;
        }}
        .results-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        .results-table th, .results-table td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #eee;
        }}
        .results-table th {{
            background: #667eea;
            color: white;
        }}
        .results-table tr:hover {{
            background: #f8f8f8;
        }}
        .disclaimer {{
            background: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 20px;
            margin: 20px 0;
            border-radius: 5px;
        }}
        .disclaimer h3 {{
            color: #856404;
            margin-bottom: 10px;
        }}
        .chart-container {{
            margin: 20px 0;
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            font-size: 0.9em;
        }}
        .download-btn {{
            display: inline-block;
            background: #667eea;
            color: white;
            padding: 10px 20px;
            border-radius: 5px;
            text-decoration: none;
            margin: 10px 5px;
        }}
        .download-btn:hover {{
            background: #5a6fd6;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 DIY Ancestry Analysis Results</h1>
        <p>Analysis completed: {timestamp}</p>
    </div>
    
    <div class="container">
        <div class="disclaimer">
            <h3>⚠️ Important Disclaimer</h3>
            <p>These results are <strong>statistical estimates</strong> based on comparison with reference 
            populations. They should not be interpreted as definitive ancestry percentages. Genetic ancestry 
            analysis has inherent limitations:</p>
            <ul style="margin-left: 20px; margin-top: 10px;">
                <li>Reference panels do not represent all genetic diversity within countries</li>
                <li>Modern country borders do not reflect historical population movements</li>
                <li>Results may vary with different reference datasets and methods</li>
                <li>These estimates reflect genetic similarity, not documented genealogy</li>
            </ul>
        </div>
        
        <div class="section">
            <h2>📊 Summary Results</h2>
            {summary_table}
        </div>
        
        <div class="section">
            <h2>🗺️ Ancestry Breakdown</h2>
            <div id="ancestry-chart" class="chart-container"></div>
        </div>
        
        <div class="section">
            <h2>📈 PCA Visualization</h2>
            <p>Your position (★) relative to European reference populations.</p>
            <div id="pca-chart" class="chart-container"></div>
        </div>
        
        <div class="section">
            <h2>🧩 ADMIXTURE Components</h2>
            <p>Your ancestral component proportions from the ADMIXTURE analysis.</p>
            <div id="admixture-chart" class="chart-container"></div>
        </div>
        
        <div class="section">
            <h2>📥 Download Results</h2>
            <p>Download your results in various formats:</p>
            <a href="#" class="download-btn" onclick="downloadCSV()">Download CSV</a>
        </div>
    </div>
    
    <div class="footer">
        <p>Generated using DIY Ancestry Analysis Pipeline</p>
        <p>PCA weight: {pca_weight}% | ADMIXTURE weight: {admixture_weight}%</p>
    </div>
    
    <script>
        // Ancestry Chart
        {ancestry_chart_js}
        
        // PCA Chart
        {pca_chart_js}
        
        // ADMIXTURE Chart
        {admixture_chart_js}
        
        // Download function
        function downloadCSV() {{
            var data = {csv_data};
            var csv = 'Country,PCA_%,ADMIXTURE_%,Ensemble_%\\n';
            data.forEach(function(row) {{
                csv += row.join(',') + '\\n';
            }});
            var blob = new Blob([csv], {{type: 'text/csv'}});
            var url = window.URL.createObjectURL(blob);
            var a = document.createElement('a');
            a.href = url;
            a.download = 'ancestry_results.csv';
            a.click();
        }}
    </script>
</body>
</html>
"""


def load_ensemble_results(results_dir: str) -> Optional[pd.DataFrame]:
    """Load ensemble results."""
    ensemble_file = os.path.join(results_dir, 'FINAL_ensemble.csv')
    
    if os.path.exists(ensemble_file):
        return pd.read_csv(ensemble_file)
    
    return None


def load_pca_data(results_dir: str) -> Optional[pd.DataFrame]:
    """Load PCA data."""
    eigenvec_file = os.path.join(results_dir, 'pca', 'eigenvec.txt')
    
    if not os.path.exists(eigenvec_file):
        return None
    
    with open(eigenvec_file, 'r') as f:
        first_line = f.readline().strip()
        has_header = 'PC1' in first_line.upper()
    
    if has_header:
        return pd.read_csv(eigenvec_file, sep=r'\s+', header=0)
    else:
        cols = ['FID', 'IID'] + [f'PC{i}' for i in range(1, 25)]
        return pd.read_csv(eigenvec_file, sep=r'\s+', header=None, names=cols)


def load_admixture_components(results_dir: str) -> Optional[pd.DataFrame]:
    """Load ADMIXTURE component breakdown."""
    probs_dir = os.path.join(results_dir, 'probabilities')
    comp_file = os.path.join(probs_dir, 'admixture_probabilities_components.csv')
    
    if os.path.exists(comp_file):
        return pd.read_csv(comp_file)
    
    return None


def create_ancestry_chart(results_df: pd.DataFrame) -> str:
    """Create Plotly ancestry bar chart."""
    if results_df is None or results_df.empty:
        return "Plotly.newPlot('ancestry-chart', [], {});"
    
    # Get top countries
    top_results = results_df.nlargest(10, 'Ensemble_%')
    
    fig = go.Figure(data=[
        go.Bar(name='PCA', x=top_results['Country'], y=top_results['PCA_%'], 
               marker_color='#667eea'),
        go.Bar(name='ADMIXTURE', x=top_results['Country'], y=top_results['ADMIXTURE_%'],
               marker_color='#764ba2'),
        go.Bar(name='Ensemble', x=top_results['Country'], y=top_results['Ensemble_%'],
               marker_color='#22c55e')
    ])
    
    fig.update_layout(
        barmode='group',
        xaxis_title='Country',
        yaxis_title='Probability (%)',
        legend=dict(orientation='h', yanchor='bottom', y=1.02),
        height=400
    )
    
    return f"Plotly.newPlot('ancestry-chart', {fig.to_json()});"


def create_pca_chart(pca_df: pd.DataFrame, metadata_df: pd.DataFrame = None) -> str:
    """Create Plotly PCA scatter plot."""
    if pca_df is None or pca_df.empty:
        return "Plotly.newPlot('pca-chart', [], {});"
    
    pca_df.columns = pca_df.columns.str.upper()
    
    fig = go.Figure()
    
    # Simple version without metadata
    fig.add_trace(go.Scatter(
        x=pca_df['PC1'],
        y=pca_df['PC2'],
        mode='markers',
        marker=dict(size=8, color='steelblue', opacity=0.6),
        name='Reference'
    ))
    
    # Highlight first sample as user
    fig.add_trace(go.Scatter(
        x=[pca_df['PC1'].iloc[0]],
        y=[pca_df['PC2'].iloc[0]],
        mode='markers',
        marker=dict(size=20, color='red', symbol='star'),
        name='You'
    ))
    
    fig.update_layout(
        xaxis_title='PC1',
        yaxis_title='PC2',
        height=500,
        showlegend=True
    )
    
    return f"Plotly.newPlot('pca-chart', {fig.to_json()});"


def create_admixture_chart(components_df: pd.DataFrame) -> str:
    """Create Plotly ADMIXTURE pie chart."""
    if components_df is None or components_df.empty:
        return "Plotly.newPlot('admixture-chart', [], {});"
    
    fig = go.Figure(data=[go.Pie(
        labels=components_df['Component'],
        values=components_df['Proportion'] * 100,
        hovertemplate='%{label}: %{value:.1f}%<br>%{customdata}<extra></extra>',
        customdata=components_df['Interpretation']
    )])
    
    fig.update_layout(height=400)
    
    return f"Plotly.newPlot('admixture-chart', {fig.to_json()});"


def create_summary_table(results_df: pd.DataFrame) -> str:
    """Create HTML summary table."""
    if results_df is None or results_df.empty:
        return "<p>No results available.</p>"
    
    top_results = results_df.nlargest(10, 'Ensemble_%')
    
    rows = ""
    for _, row in top_results.iterrows():
        rows += f"""
        <tr>
            <td>{row['Country']}</td>
            <td>{row['PCA_%']:.1f}%</td>
            <td>{row['ADMIXTURE_%']:.1f}%</td>
            <td><strong>{row['Ensemble_%']:.1f}%</strong></td>
        </tr>
        """
    
    return f"""
    <table class="results-table">
        <thead>
            <tr>
                <th>Country</th>
                <th>PCA Method</th>
                <th>ADMIXTURE Method</th>
                <th>Ensemble</th>
            </tr>
        </thead>
        <tbody>
            {rows}
        </tbody>
    </table>
    """


def create_csv_data(results_df: pd.DataFrame) -> str:
    """Create CSV data for download."""
    if results_df is None or results_df.empty:
        return "[]"
    
    data = []
    for _, row in results_df.iterrows():
        data.append([
            f'"{row["Country"]}"',
            f'{row["PCA_%"]:.2f}',
            f'{row["ADMIXTURE_%"]:.2f}',
            f'{row["Ensemble_%"]:.2f}'
        ])
    
    return str(data)


def generate_report(results_dir: str, output_file: str):
    """Generate the complete HTML report."""
    logger.info("Generating HTML report...")
    
    # Load data
    results_df = load_ensemble_results(results_dir)
    pca_df = load_pca_data(results_dir)
    components_df = load_admixture_components(results_dir)
    
    # Create charts
    ancestry_chart_js = create_ancestry_chart(results_df)
    pca_chart_js = create_pca_chart(pca_df)
    admixture_chart_js = create_admixture_chart(components_df)
    
    # Create table and data
    summary_table = create_summary_table(results_df)
    csv_data = create_csv_data(results_df)
    
    # Generate HTML
    html_content = HTML_TEMPLATE.format(
        timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        summary_table=summary_table,
        ancestry_chart_js=ancestry_chart_js,
        pca_chart_js=pca_chart_js,
        admixture_chart_js=admixture_chart_js,
        csv_data=csv_data,
        pca_weight=30,
        admixture_weight=70
    )
    
    # Write output
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    logger.info(f"Report saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate interactive HTML dashboard',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--results-dir', '-r',
        required=True,
        help='Path to results directory'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='interactive_dashboard.html',
        help='Output HTML file'
    )
    
    args = parser.parse_args()
    
    generate_report(args.results_dir, args.output)


if __name__ == '__main__':
    main()
