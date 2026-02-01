"""
Country Predictor Module
Combines PCA and ADMIXTURE results to predict country-level ancestry with uncertainty.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from collections import defaultdict


class CountryPredictor:
    """
    Predicts country-level ancestry percentages with uncertainty estimates.
    Combines results from PCA and ADMIXTURE analyses.
    """
    
    def __init__(self, uncertainty_threshold: float = 5.0):
        """
        Initialize country predictor.
        
        Args:
            uncertainty_threshold: Minimum percentage to report (filters noise)
        """
        self.uncertainty_threshold = uncertainty_threshold
        
    def combine_predictions(self, pca_results: Dict[str, float],
                          admixture_results: Dict[str, float],
                          pca_weight: float = 0.5) -> Dict[str, float]:
        """
        Combine PCA and ADMIXTURE predictions into final ancestry estimates.
        
        Args:
            pca_results: PCA-based ancestry percentages
            admixture_results: ADMIXTURE-based ancestry percentages
            pca_weight: Weight for PCA results (0-1), rest goes to ADMIXTURE
            
        Returns:
            Combined ancestry percentages
        """
        combined = {}
        all_populations = set(pca_results.keys()) | set(admixture_results.keys())
        
        admixture_weight = 1.0 - pca_weight
        
        for pop in all_populations:
            pca_val = pca_results.get(pop, 0.0)
            admix_val = admixture_results.get(pop, 0.0)
            
            combined[pop] = (pca_val * pca_weight) + (admix_val * admixture_weight)
        
        # Normalize to 100%
        total = sum(combined.values())
        if total > 0:
            combined = {pop: (val / total) * 100 for pop, val in combined.items()}
        
        # Filter out small contributions
        combined = {pop: val for pop, val in combined.items() 
                   if val >= self.uncertainty_threshold}
        
        # Re-normalize after filtering
        total = sum(combined.values())
        if total > 0:
            combined = {pop: (val / total) * 100 for pop, val in combined.items()}
        
        return combined
    
    def calculate_uncertainty(self, pca_results: Dict[str, float],
                            admixture_results: Dict[str, float]) -> Dict[str, Dict[str, float]]:
        """
        Calculate uncertainty metrics for each population prediction.
        
        Args:
            pca_results: PCA-based ancestry percentages
            admixture_results: ADMIXTURE-based ancestry percentages
            
        Returns:
            Dictionary with uncertainty metrics per population
        """
        uncertainty = {}
        all_populations = set(pca_results.keys()) | set(admixture_results.keys())
        
        for pop in all_populations:
            pca_val = pca_results.get(pop, 0.0)
            admix_val = admixture_results.get(pop, 0.0)
            
            # Calculate disagreement between methods
            disagreement = abs(pca_val - admix_val)
            
            # Calculate average and range
            values = [v for v in [pca_val, admix_val] if v > 0]
            if values:
                avg = np.mean(values)
                range_val = max(values) - min(values)
                
                uncertainty[pop] = {
                    'mean': avg,
                    'disagreement': disagreement,
                    'range': range_val,
                    'pca_estimate': pca_val,
                    'admixture_estimate': admix_val
                }
        
        return uncertainty
    
    def identify_regional_clusters(self, predictions: Dict[str, float]) -> Dict[str, List[str]]:
        """
        Group predicted ancestries into geographic regions.
        
        Args:
            predictions: Country-level ancestry predictions
            
        Returns:
            Dictionary mapping regions to list of countries
        """
        # Define European regions
        regions = {
            'Western Europe': ['British', 'Irish', 'French', 'Dutch', 'Belgian'],
            'Southern Europe': ['Spanish', 'Italian', 'Portuguese', 'Greek'],
            'Northern Europe': ['Swedish', 'Norwegian', 'Danish', 'Finnish', 'Estonian', 'Latvian', 'Lithuanian'],
            'Eastern Europe': ['Polish', 'Romanian', 'Czech', 'Hungarian', 'Bulgarian', 
                              'Croatian', 'Serbian', 'Slovenian', 'Slovak', 'Ukrainian', 'Russian'],
            'Central Europe': ['German', 'Austrian', 'Swiss']
        }
        
        regional_breakdown = defaultdict(list)
        
        for country, percentage in predictions.items():
            # Find which region this country belongs to
            for region, countries in regions.items():
                if country in countries:
                    regional_breakdown[region].append((country, percentage))
                    break
        
        # Sort countries within each region by percentage
        for region in regional_breakdown:
            regional_breakdown[region].sort(key=lambda x: x[1], reverse=True)
        
        return dict(regional_breakdown)
    
    def generate_report(self, combined_predictions: Dict[str, float],
                       uncertainty: Dict[str, Dict[str, float]],
                       confidence_metrics: Optional[Dict[str, float]] = None) -> str:
        """
        Generate a human-readable ancestry report.
        
        Args:
            combined_predictions: Final ancestry predictions
            uncertainty: Uncertainty metrics
            confidence_metrics: Optional confidence scores
            
        Returns:
            Formatted report string
        """
        report_lines = []
        report_lines.append("=" * 60)
        report_lines.append("EUROPEAN ANCESTRY ANALYSIS REPORT")
        report_lines.append("=" * 60)
        report_lines.append("")
        
        # Add confidence warning
        report_lines.append("NOTE: This is a DIY ancestry analysis tool.")
        report_lines.append("Results are estimates and should be interpreted with caution.")
        report_lines.append("Uncertainty and overlap are especially high in regions like")
        report_lines.append("the Balkans and Central Europe.")
        report_lines.append("")
        
        # Sort predictions by percentage
        sorted_preds = sorted(combined_predictions.items(), 
                            key=lambda x: x[1], reverse=True)
        
        report_lines.append("ANCESTRY PREDICTIONS:")
        report_lines.append("-" * 60)
        
        for pop, pct in sorted_preds:
            # Get uncertainty info
            unc = uncertainty.get(pop, {})
            disagreement = unc.get('disagreement', 0)
            range_val = unc.get('range', 0)
            
            report_lines.append(f"{pop:20s}: {pct:5.1f}%")
            
            # Add uncertainty indicator
            if disagreement > 10:
                report_lines.append(f"{'':20s}  ⚠ High uncertainty (±{disagreement:.1f}%)")
            elif range_val > 5:
                report_lines.append(f"{'':20s}  ⚠ Moderate uncertainty (±{range_val:.1f}%)")
            
            # Show method-specific estimates if they differ significantly
            if disagreement > 10:
                pca_est = unc.get('pca_estimate', 0)
                admix_est = unc.get('admixture_estimate', 0)
                report_lines.append(
                    f"{'':20s}    PCA: {pca_est:.1f}%, ADMIXTURE: {admix_est:.1f}%"
                )
        
        report_lines.append("")
        
        # Regional summary
        regional = self.identify_regional_clusters(combined_predictions)
        if regional:
            report_lines.append("REGIONAL SUMMARY:")
            report_lines.append("-" * 60)
            
            for region, countries in regional.items():
                total_pct = sum(pct for _, pct in countries)
                report_lines.append(f"{region:30s}: {total_pct:5.1f}%")
                for country, pct in countries:
                    report_lines.append(f"  - {country:25s}: {pct:5.1f}%")
            
            report_lines.append("")
        
        # Add overall confidence if available
        if confidence_metrics:
            report_lines.append("CONFIDENCE METRICS:")
            report_lines.append("-" * 60)
            overall_conf = confidence_metrics.get('overall_confidence', 0)
            report_lines.append(f"Overall Confidence: {overall_conf:.2f} (0-1 scale)")
            
            if overall_conf < 0.5:
                report_lines.append("⚠ Low confidence - results may be highly uncertain")
            elif overall_conf < 0.7:
                report_lines.append("⚠ Moderate confidence - interpret with caution")
            else:
                report_lines.append("✓ Good confidence in results")
        
        report_lines.append("")
        report_lines.append("=" * 60)
        report_lines.append("For more accurate results, consider professional testing")
        report_lines.append("or using larger, more diverse reference datasets.")
        report_lines.append("=" * 60)
        
        return "\n".join(report_lines)
    
    def get_top_matches(self, predictions: Dict[str, float], 
                       n: int = 5) -> List[Tuple[str, float]]:
        """
        Get top N ancestry matches.
        
        Args:
            predictions: Ancestry predictions
            n: Number of top matches to return
            
        Returns:
            List of (country, percentage) tuples
        """
        sorted_preds = sorted(predictions.items(), key=lambda x: x[1], reverse=True)
        return sorted_preds[:n]
