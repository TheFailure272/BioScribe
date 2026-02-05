import requests
import json

# Test the backend API
try:
    print("Testing backend connection...")
    
    # Test health endpoint
    health_response = requests.get("http://localhost:8000/api/health")
    print(f"Health check: {health_response.status_code}")
    print(f"Health response: {health_response.json()}")
    
    # Test protein analysis
    test_data = {
        "sequence": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "name": "Test Protein",
        "organism": "Human"
    }
    
    print("\nTesting protein analysis...")
    analysis_response = requests.post(
        "http://localhost:8000/api/ai/analyze-protein",
        json=test_data,
        headers={"Content-Type": "application/json"}
    )
    
    print(f"Analysis status: {analysis_response.status_code}")
    if analysis_response.status_code == 200:
        result = analysis_response.json()
        print(f"Analysis successful!")
        print(f"Protein name: {result.get('name')}")
        print(f"Sequence length: {result.get('length')}")
        print(f"Molecular weight: {result.get('molecular_properties', {}).get('molecular_weight')}")
    else:
        print(f"Analysis failed: {analysis_response.text}")
        
except Exception as e:
    print(f"Error testing backend: {e}")
