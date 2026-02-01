"""
PCA Analysis Module
Projects genomes into European genetic space using Principal Component Analysis.
"""

import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from typing import Tuple, Optional, Dict
import warnings


class PCAAnalyzer:
    """
    Performs PCA on genetic data to project individuals into genetic space.
    """
    
    def __init__(self, n_components: int = 10):
        """
        Initialize PCA analyzer.
        
        Args:
            n_components: Number of principal components to compute
        """
        self.n_components = n_components
        self.pca = None
        self.scaler = StandardScaler()
        self.reference_projections = None
        self.reference_labels = None
        
    def fit_reference(self, reference_data: np.ndarray, 
                     reference_labels: Optional[np.ndarray] = None):
        """
        Fit PCA model on reference data.
        
        Args:
            reference_data: Reference genetic data matrix (samples x SNPs)
            reference_labels: Population labels for reference samples
        """
        # Standardize the data
        scaled_data = self.scaler.fit_transform(reference_data)
        
        # Fit PCA
        self.pca = PCA(n_components=self.n_components)
        self.reference_projections = self.pca.fit_transform(scaled_data)
        self.reference_labels = reference_labels
        
    def transform_user(self, user_data: np.ndarray) -> np.ndarray:
        """
        Project user data into the reference PCA space.
        
        Args:
            user_data: User genetic data matrix (samples x SNPs)
            
        Returns:
            Projected coordinates in PCA space
        """
        if self.pca is None:
            raise ValueError("PCA model not fitted. Call fit_reference first.")
        
        # Ensure user data has same number of features
        if user_data.shape[1] != self.pca.n_features_in_:
            raise ValueError(
                f"User data has {user_data.shape[1]} features, "
                f"but model expects {self.pca.n_features_in_}"
            )
        
        scaled_data = self.scaler.transform(user_data)
        user_projection = self.pca.transform(scaled_data)
        
        return user_projection
    
    def get_explained_variance(self) -> np.ndarray:
        """
        Get the explained variance ratio for each component.
        
        Returns:
            Array of explained variance ratios
        """
        if self.pca is None:
            raise ValueError("PCA model not fitted.")
        
        return self.pca.explained_variance_ratio_
    
    def get_nearest_populations(self, user_projection: np.ndarray, 
                               k: int = 5) -> Dict[str, float]:
        """
        Find k nearest reference populations to user in PCA space.
        
        Args:
            user_projection: User's coordinates in PCA space
            k: Number of nearest neighbors to find
            
        Returns:
            Dictionary mapping population names to distances
        """
        if self.reference_projections is None or self.reference_labels is None:
            raise ValueError("Reference data not fitted.")
        
        # Calculate distances to all reference samples
        distances = np.sqrt(np.sum((self.reference_projections - user_projection) ** 2, axis=1))
        
        # Find k nearest
        nearest_indices = np.argsort(distances)[:k]
        nearest_labels = self.reference_labels[nearest_indices]
        nearest_distances = distances[nearest_indices]
        
        # Group by population and take mean distance
        pop_distances = {}
        for label, dist in zip(nearest_labels, nearest_distances):
            if label not in pop_distances:
                pop_distances[label] = []
            pop_distances[label].append(dist)
        
        # Average distances per population
        pop_distances = {pop: np.mean(dists) for pop, dists in pop_distances.items()}
        
        return pop_distances
    
    def calculate_population_proximity(self, user_projection: np.ndarray, 
                                      n_neighbors: int = 50) -> Dict[str, float]:
        """
        Calculate proximity to each population based on PCA space.
        
        Args:
            user_projection: User's coordinates in PCA space
            n_neighbors: Number of neighbors to consider
            
        Returns:
            Dictionary mapping population names to proximity scores (0-1)
        """
        if self.reference_projections is None or self.reference_labels is None:
            raise ValueError("Reference data not fitted.")
        
        # Calculate distances to all reference samples
        distances = np.sqrt(np.sum((self.reference_projections - user_projection) ** 2, axis=1))
        
        # Find n nearest neighbors
        nearest_indices = np.argsort(distances)[:n_neighbors]
        nearest_labels = self.reference_labels[nearest_indices]
        nearest_distances = distances[nearest_indices]
        
        # Convert distances to similarities (inverse of distance)
        # Use exponential decay for more natural weighting
        similarities = np.exp(-nearest_distances / np.mean(nearest_distances))
        
        # Count weighted votes for each population
        pop_scores = {}
        for label, similarity in zip(nearest_labels, similarities):
            pop_scores[label] = pop_scores.get(label, 0) + similarity
        
        # Normalize to percentages
        total_score = sum(pop_scores.values())
        pop_percentages = {pop: (score / total_score) * 100 
                          for pop, score in pop_scores.items()}
        
        return pop_percentages
    
    def get_component_loadings(self, n_top_snps: int = 10) -> Dict[int, list]:
        """
        Get top SNP loadings for each principal component.
        
        Args:
            n_top_snps: Number of top SNPs to return per component
            
        Returns:
            Dictionary mapping component index to list of (SNP_index, loading) tuples
        """
        if self.pca is None:
            raise ValueError("PCA model not fitted.")
        
        loadings = {}
        for i in range(self.n_components):
            component = self.pca.components_[i]
            top_indices = np.argsort(np.abs(component))[-n_top_snps:][::-1]
            loadings[i] = [(idx, component[idx]) for idx in top_indices]
        
        return loadings
