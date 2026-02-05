// BioScribe AI API Client
const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

export interface ProteinAnalysisRequest {
  sequence: string;
  name?: string;
  organism?: string;
}

export interface DrugGenerationRequest {
  protein_id: string;
  num_candidates?: number;
}

export interface DockingRequest {
  session_id: string;
}

// New AI-Enhanced Interfaces
export interface AIProteinAnalysisRequest {
  sequence: string;
  name?: string;
  organism?: string;
}

export interface AIMoleculeGenerationRequest {
  sequence: string;
  name?: string;
  organism?: string;
  num_molecules?: number;
  target_properties?: Record<string, number>;
}

export interface AIEnhancedDockingRequest {
  session_id: string;
  num_poses?: number;
  use_alphafold?: boolean;
}

export interface AICompletePipelineRequest {
  sequence: string;
  name?: string;
  organism?: string;
  num_molecules?: number;
  num_poses?: number;
  target_properties?: Record<string, number>;
}

class BioScribeAPI {
  private baseUrl: string;

  constructor(baseUrl: string = API_BASE_URL) {
    this.baseUrl = baseUrl;
  }

  private async request<T>(
    endpoint: string,
    options: RequestInit = {}
  ): Promise<T> {
    const url = `${this.baseUrl}${endpoint}`;
    
    const config: RequestInit = {
      headers: {
        'Content-Type': 'application/json',
        ...options.headers,
      },
      ...options,
    };

    try {
      const response = await fetch(url, config);
      
      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        throw new Error(errorData.detail || `HTTP error! status: ${response.status}`);
      }

      return await response.json();
    } catch (error) {
      console.error(`API request failed: ${endpoint}`, error);
      throw error;
    }
  }

  // Health check
  async healthCheck() {
    return this.request('/health');
  }

  // Get API status
  async getStatus() {
    return this.request('/');
  }

  // Protein analysis
  async analyzeProtein(data: ProteinAnalysisRequest) {
    return this.request('/api/protein/analyze', {
      method: 'POST',
      body: JSON.stringify(data),
    });
  }

  // Drug generation
  async generateDrugs(data: DrugGenerationRequest) {
    return this.request('/api/drugs/generate', {
      method: 'POST',
      body: JSON.stringify(data),
    });
  }

  // Docking simulation
  async startDocking(data: DockingRequest) {
    return this.request('/api/docking/simulate', {
      method: 'POST',
      body: JSON.stringify(data),
    });
  }

  // Get docking status
  async getDockingStatus(sessionId: string) {
    return this.request(`/api/docking/status/${sessionId}`);
  }
  // Get results
  async getResults(sessionId: string) {
    return this.request(`/api/results/${sessionId}`);
  }

  // Export methods
  async exportReport(sessionId: string, format: 'json' | 'report' | 'summary' | 'text') {
    return this.request(`/api/export/${sessionId}/${format}`);
  }

  // AI-Enhanced Methods
  async aiAnalyzeProtein(request: AIProteinAnalysisRequest) {
    return this.request('/api/ai/analyze-protein', {
      method: 'POST',
      body: JSON.stringify(request),
    });
  }

  async aiGenerateMolecules(request: AIMoleculeGenerationRequest) {
    return this.request('/api/ai/generate-molecules', {
      method: 'POST',
      body: JSON.stringify(request),
    });
  }

  async aiEnhancedDocking(request: AIEnhancedDockingRequest) {
    return this.request('/api/ai/enhanced-docking', {
      method: 'POST',
      body: JSON.stringify(request),
    });
  }

  async aiCompletePipeline(request: AICompletePipelineRequest) {
    return this.request('/api/ai/complete-pipeline', {
      method: 'POST',
      body: JSON.stringify(request),
    });
  }

  async getAISession(sessionId: string) {
    return this.request(`/api/ai/session/${sessionId}`);
  }

  async checkAIHealth() {
    return this.request('/api/ai/health');
  }

  // Export comprehensive medical report
  async exportMedicalReport(sessionId: string) {
    return this.request(`/api/export/${sessionId}/report`);
  }
  async exportSummary(sessionId: string) {
    return this.request(`/api/export/${sessionId}/summary`);
  }

  // Export as text format
  async exportTextReport(sessionId: string) {
    return this.request(`/api/export/${sessionId}/text`);
  }
}

// Create singleton instance
export const apiClient = new BioScribeAPI();

// Helper functions for common operations
export const BioScribeHelpers = {
  // Check if API is available
  async isApiAvailable(): Promise<boolean> {
    try {
      await apiClient.healthCheck();
      return true;
    } catch {
      return false;
    }
  },

  // Format API errors for display
  formatError(error: unknown): string {
    if (error instanceof Error) {
      return error.message;
    }
    return 'An unexpected error occurred';
  },

  // Validate protein sequence
  validateProteinSequence(sequence: string): { isValid: boolean; error?: string } {
    const cleanSeq = sequence.replace(/\s+/g, '').toUpperCase();
    
    if (cleanSeq.length < 10) {
      return { isValid: false, error: 'Sequence too short (minimum 10 amino acids)' };
    }
    
    const validAminoAcids = /^[ACDEFGHIKLMNPQRSTVWY]+$/;
    if (!validAminoAcids.test(cleanSeq)) {
      return { isValid: false, error: 'Invalid amino acid characters found' };
    }
    
    return { isValid: true };
  },

  // Generate mock data when API is unavailable
  generateMockProteinAnalysis: (sequence: string, name?: string) => ({
    protein_id: `mock_${Date.now()}`,
    name: name || 'Mock Protein',
    sequence: sequence.replace(/\s+/g, '').toUpperCase(),
    length: sequence.replace(/\s+/g, '').length,
    molecular_weight: 25000 + Math.random() * 50000,
    properties: {
      isoelectric_point: 6.0 + Math.random() * 4,
      instability_index: 20 + Math.random() * 40,
      gravy: -0.5 + Math.random(),
      aromaticity: Math.random() * 0.2,
    },
    druggability_score: 0.7 + Math.random() * 0.3,
    analysis_timestamp: new Date().toISOString(),
  }),

  generateMockDrugCandidates: (count: number = 10) => {
    const candidates = [];
    const templates = [
      'CC1=CC(=O)NC(=O)N1',
      'COc1cc2c(cc1OC)c(=O)c(cn2)c3cccc(c3)Cl',
      'Cc1ccc(cc1)c2cc(nn2c3ccc(cc3)S(=O)(=O)N)C(F)(F)F',
    ];

    for (let i = 0; i < count; i++) {
      candidates.push({
        smiles: templates[i % templates.length],
        name: `MockDrug-${String(i + 1).padStart(3, '0')}`,
        molecular_weight: 200 + Math.random() * 400,
        logP: -1 + Math.random() * 6,
        tpsa: 20 + Math.random() * 120,
        qed: 0.3 + Math.random() * 0.7,
        drug_likeness_score: 0.5 + Math.random() * 0.5,
      });
    }

    return candidates.sort((a, b) => b.drug_likeness_score - a.drug_likeness_score);
  },

  generateMockDockingResults: (candidates: any[]) => {
    return candidates.map((candidate, index) => ({
      ...candidate,
      binding_affinity: -12 + Math.random() * 8,
      rmsd: 0.5 + Math.random() * 3,
      confidence: 60 + Math.random() * 35,
      interactions: [
        { type: 'H-bond', residue: 'Asp25', distance: 2.1 + Math.random() * 0.8 },
        { type: 'Hydrophobic', residue: 'Ile50', distance: 3.2 + Math.random() * 1.0 },
      ].slice(0, Math.floor(1 + Math.random() * 2))
    })).sort((a, b) => a.binding_affinity - b.binding_affinity);
  }
};

export default apiClient;
export const api = apiClient;
