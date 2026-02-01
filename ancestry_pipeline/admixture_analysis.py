"""
ADMIXTURE Analysis Module
Estimates ancestry proportions using ADMIXTURE-like algorithms.
"""

import numpy as np
from sklearn.cluster import KMeans
from typing import Dict, List, Optional, Tuple
import warnings


class AdmixtureAnalyzer:
    """
    Performs ADMIXTURE-like analysis to estimate ancestry proportions.
    Uses a simplified approach based on clustering and proximity.
    """
    
    def __init__(self, n_clusters: int = 5, max_iterations: int = 100):
        """
        Initialize ADMIXTURE analyzer.
        
        Args:
            n_clusters: Number of ancestral populations (K)
            max_iterations: Maximum iterations for optimization
        """
        self.n_clusters = n_clusters
        self.max_iterations = max_iterations
        self.cluster_model = None
        self.population_clusters = None
        self.reference_labels = None
        
    def fit_reference(self, reference_data: np.ndarray, 
                     reference_labels: Optional[np.ndarray] = None):
        """
        Fit ADMIXTURE model on reference data.
        
        Args:
            reference_data: Reference genetic data matrix (samples x SNPs)
            reference_labels: Population labels for reference samples
        """
        self.reference_labels = reference_labels
        
        # Use K-means clustering as a simplified approach
        self.cluster_model = KMeans(n_clusters=self.n_clusters, 
                                   max_iter=self.max_iterations,
                                   n_init=10,
                                   random_state=42)
        
        cluster_assignments = self.cluster_model.fit_predict(reference_data)
        
        # Map clusters to populations if labels provided
        if reference_labels is not None:
            self.population_clusters = self._map_clusters_to_populations(
                cluster_assignments, reference_labels
            )
        
    def _map_clusters_to_populations(self, cluster_assignments: np.ndarray,
                                     labels: np.ndarray) -> Dict[int, Dict[str, float]]:
        """
        Map clusters to population distributions.
        
        Args:
            cluster_assignments: Cluster assignment for each sample
            labels: Population labels
            
        Returns:
            Dictionary mapping cluster ID to population distribution
        """
        cluster_pops = {}
        
        for cluster_id in range(self.n_clusters):
            cluster_mask = cluster_assignments == cluster_id
            cluster_labels = labels[cluster_mask]
            
            if len(cluster_labels) > 0:
                unique_pops, counts = np.unique(cluster_labels, return_counts=True)
                total = counts.sum()
                pop_dist = {pop: (count / total) * 100 
                           for pop, count in zip(unique_pops, counts)}
                cluster_pops[cluster_id] = pop_dist
            else:
                cluster_pops[cluster_id] = {}
        
        return cluster_pops
    
    def estimate_ancestry(self, user_data: np.ndarray) -> Dict[str, float]:
        """
        Estimate ancestry proportions for user data.
        
        Args:
            user_data: User genetic data matrix (samples x SNPs)
            
        Returns:
            Dictionary mapping population names to ancestry percentages
        """
        if self.cluster_model is None:
            raise ValueError("Model not fitted. Call fit_reference first.")
        
        # Get distances to each cluster centroid
        distances = self.cluster_model.transform(user_data)
        
        # Convert distances to probabilities using softmax-like approach
        # Invert distances so closer = higher weight
        max_dist = distances.max()
        inv_distances = max_dist - distances
        
        # Apply exponential weighting
        weights = np.exp(inv_distances / np.mean(inv_distances))
        
        # Normalize to get cluster membership probabilities
        cluster_probs = weights / weights.sum(axis=1, keepdims=True)
        
        # If we have population mappings, convert cluster probs to population probs
        if self.population_clusters is not None:
            ancestry = self._cluster_probs_to_ancestry(cluster_probs[0])
        else:
            # Without population labels, return cluster memberships
            ancestry = {f'Cluster_{i}': prob * 100 
                       for i, prob in enumerate(cluster_probs[0])}
        
        return ancestry
    
    def _cluster_probs_to_ancestry(self, cluster_probs: np.ndarray) -> Dict[str, float]:
        """
        Convert cluster membership probabilities to ancestry percentages.
        
        Args:
            cluster_probs: Probability of belonging to each cluster
            
        Returns:
            Dictionary mapping population names to ancestry percentages
        """
        ancestry = {}
        
        for cluster_id, cluster_prob in enumerate(cluster_probs):
            if cluster_id in self.population_clusters:
                pop_dist = self.population_clusters[cluster_id]
                
                # Add weighted contribution of this cluster to each population
                for pop, pop_pct_in_cluster in pop_dist.items():
                    contribution = cluster_prob * (pop_pct_in_cluster / 100)
                    ancestry[pop] = ancestry.get(pop, 0) + contribution
        
        # Convert to percentages
        total = sum(ancestry.values())
        if total > 0:
            ancestry = {pop: (val / total) * 100 for pop, val in ancestry.items()}
        
        return ancestry
    
    def get_cluster_info(self) -> Dict[int, Dict]:
        """
        Get information about each cluster.
        
        Returns:
            Dictionary with cluster statistics and population distributions
        """
        if self.cluster_model is None:
            raise ValueError("Model not fitted.")
        
        cluster_info = {}
        for cluster_id in range(self.n_clusters):
            info = {
                'centroid': self.cluster_model.cluster_centers_[cluster_id],
                'inertia': self.cluster_model.inertia_,
            }
            
            if self.population_clusters and cluster_id in self.population_clusters:
                info['population_distribution'] = self.population_clusters[cluster_id]
            
            cluster_info[cluster_id] = info
        
        return cluster_info
    
    def estimate_admixture_confidence(self, user_data: np.ndarray) -> Dict[str, float]:
        """
        Estimate confidence scores for ancestry estimates.
        
        Args:
            user_data: User genetic data matrix
            
        Returns:
            Dictionary with confidence metrics
        """
        if self.cluster_model is None:
            raise ValueError("Model not fitted.")
        
        distances = self.cluster_model.transform(user_data)
        
        # Calculate confidence metrics
        min_dist = distances.min()
        max_dist = distances.max()
        mean_dist = distances.mean()
        
        # Distance to nearest cluster (lower = more confident)
        nearest_cluster_dist = min_dist
        
        # Spread of distances (higher spread = more uncertain)
        dist_std = distances.std()
        
        # Simple confidence score (0-1, higher = more confident)
        # Based on how clearly user fits into one cluster vs being intermediate
        confidence = 1.0 / (1.0 + dist_std / mean_dist)
        
        return {
            'overall_confidence': confidence,
            'nearest_cluster_distance': nearest_cluster_dist,
            'distance_spread': dist_std,
            'mean_distance': mean_dist
        }
