#!/usr/bin/env python3
"""
Setup script for the Ancestry Analysis Pipeline.
Verifies installation and runs basic tests.
"""

import sys
import subprocess
import importlib


def check_python_version():
    """Check if Python version is adequate."""
    print("Checking Python version...")
    version = sys.version_info
    if version.major < 3 or (version.major == 3 and version.minor < 7):
        print(f"✗ Python 3.7+ required, but found {version.major}.{version.minor}")
        return False
    print(f"✓ Python {version.major}.{version.minor}.{version.micro}")
    return True


def install_dependencies():
    """Install required packages."""
    print("\nInstalling dependencies...")
    try:
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", "-q", "-r", "requirements.txt"
        ])
        print("✓ Dependencies installed")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ Failed to install dependencies: {e}")
        return False


def verify_imports():
    """Verify all required packages can be imported."""
    print("\nVerifying package imports...")
    required_packages = [
        'numpy',
        'pandas',
        'sklearn',
        'matplotlib',
        'scipy'
    ]
    
    all_ok = True
    for package in required_packages:
        try:
            importlib.import_module(package)
            print(f"✓ {package}")
        except ImportError:
            print(f"✗ {package} not found")
            all_ok = False
    
    return all_ok


def verify_pipeline_modules():
    """Verify pipeline modules can be imported."""
    print("\nVerifying pipeline modules...")
    
    # Add current directory to path
    sys.path.insert(0, '.')
    
    modules = [
        'ancestry_pipeline',
        'ancestry_pipeline.genotype_parser',
        'ancestry_pipeline.reference_data',
        'ancestry_pipeline.pca_analysis',
        'ancestry_pipeline.admixture_analysis',
        'ancestry_pipeline.country_predictor'
    ]
    
    all_ok = True
    for module in modules:
        try:
            importlib.import_module(module)
            print(f"✓ {module}")
        except ImportError as e:
            print(f"✗ {module}: {e}")
            all_ok = False
    
    return all_ok


def run_basic_test():
    """Run a basic test of the pipeline."""
    print("\nRunning basic pipeline test...")
    
    try:
        sys.path.insert(0, '.')
        from ancestry_pipeline import ReferenceDataHandler
        import numpy as np
        
        # Test synthetic data generation
        ref_handler = ReferenceDataHandler()
        ref_data, ref_labels = ref_handler.generate_synthetic_reference(
            n_samples_per_pop=5,
            n_snps=100
        )
        
        assert ref_data.shape[0] > 0, "No reference data generated"
        assert ref_data.shape[1] == 100, "Wrong number of SNPs"
        assert len(ref_labels) == ref_data.shape[0], "Label count mismatch"
        
        print("✓ Synthetic data generation works")
        
        # Test genotype parser with sample file
        from ancestry_pipeline import GenotypeParser
        try:
            parser = GenotypeParser()
            data = parser.parse_file('examples/sample_genotype.txt')
            print(f"✓ Genotype parsing works ({len(data)} SNPs)")
        except FileNotFoundError:
            print("⚠ Sample genotype file not found (optional)")
        
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main setup routine."""
    print("=" * 60)
    print("ANCESTRY ANALYSIS PIPELINE - SETUP")
    print("=" * 60)
    
    steps = [
        ("Python Version Check", check_python_version),
        ("Install Dependencies", install_dependencies),
        ("Verify Imports", verify_imports),
        ("Verify Pipeline Modules", verify_pipeline_modules),
        ("Run Basic Test", run_basic_test)
    ]
    
    for step_name, step_func in steps:
        if not step_func():
            print(f"\n✗ Setup failed at: {step_name}")
            print("\nPlease fix the errors above and run setup.py again.")
            sys.exit(1)
    
    print("\n" + "=" * 60)
    print("✓ SETUP COMPLETE!")
    print("=" * 60)
    print("\nYou can now use the ancestry analysis pipeline:")
    print("  python analyze_ancestry.py examples/sample_genotype.txt")
    print("\nFor more information, see README.md and USAGE.md")
    print("=" * 60)


if __name__ == '__main__':
    main()
