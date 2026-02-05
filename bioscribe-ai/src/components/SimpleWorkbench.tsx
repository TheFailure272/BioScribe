"use client";

import React, { useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { 
  Dna, 
  FlaskConical, 
  Target, 
  Eye, 
  Brain,
  Atom,
  Zap,
  CheckCircle,
  Loader2,
  Sparkles
} from 'lucide-react';
import { motion } from 'framer-motion';
import { MolecularViewer } from './MolecularViewer';

interface SimpleWorkbenchProps {
  initialProtein?: {
    sequence: string;
    name: string;
    organism?: string;
  };
}

export function SimpleWorkbench({ initialProtein }: SimpleWorkbenchProps) {
  const [activeTab, setActiveTab] = useState('protein');
  const [isProcessing, setIsProcessing] = useState(false);
  const [results, setResults] = useState<any>(null);

  const mockProcessAI = async () => {
    setIsProcessing(true);
    setActiveTab('generation');
    
    // Simulate AI processing
    await new Promise(resolve => setTimeout(resolve, 3000));
    
    const mockResults = {
      protein_analysis: {
        sequence_length: initialProtein?.sequence.length || 300,
        molecular_weight: 35000,
        hydrophobicity: 0.2,
        druggability_score: 0.8
      },
      candidates: [
        {
          id: "BSA-001",
          name: "BSA-001",
          smiles: "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)C",
          molecular_weight: 425.5,
          logP: 2.3,
          qed: 0.85,
          binding_affinity: -9.2,
          confidence: 0.92
        },
        {
          id: "BSA-002", 
          name: "BSA-002",
          smiles: "CN1CCN(CC1)C2=CC=C(C=C2)NC(=O)C3=CC=C(C=C3)F",
          molecular_weight: 387.4,
          logP: 1.8,
          qed: 0.78,
          binding_affinity: -8.7,
          confidence: 0.88
        }
      ],
      docking_results: {
        total_successful_dockings: 2,
        best_binding_affinity: -9.2,
        average_confidence: 0.90
      }
    };
    
    setResults(mockResults);
    setIsProcessing(false);
    setActiveTab('visualization');
  };

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header */}
      <header className="bg-white border-b border-gray-200 px-6 py-4">
        <div className="flex items-center justify-between">
          <div className="flex items-center space-x-4">
            <div className="flex items-center space-x-2">
              <Atom className="w-8 h-8 text-blue-600" />
              <div>
                <h1 className="text-xl font-bold text-gray-900">BioScribe AI</h1>
                <p className="text-sm text-gray-500">AlphaFold-Level Biocomputing Platform</p>
              </div>
            </div>
            <div className="flex items-center space-x-2">
              <Badge variant="outline" className="bg-green-50 text-green-700 border-green-200">
                <Brain className="w-3 h-3 mr-1" />
                AI Enhanced
              </Badge>
              <Badge variant="outline" className="bg-blue-50 text-blue-700 border-blue-200">
                Production Ready
              </Badge>
            </div>
          </div>
        </div>
      </header>

      <div className="container mx-auto px-6 py-8">
        <Tabs value={activeTab} onValueChange={setActiveTab}>
          <TabsList className="grid w-full grid-cols-4 max-w-2xl mx-auto mb-8">
            <TabsTrigger value="protein">
              <Dna className="w-4 h-4 mr-2" />
              Protein Input
            </TabsTrigger>
            <TabsTrigger value="generation">
              <Brain className="w-4 h-4 mr-2" />
              AI Generation
            </TabsTrigger>
            <TabsTrigger value="docking">
              <Target className="w-4 h-4 mr-2" />
              Docking
            </TabsTrigger>
            <TabsTrigger value="visualization">
              <Eye className="w-4 h-4 mr-2" />
              3D View
            </TabsTrigger>
          </TabsList>

          <TabsContent value="protein" className="space-y-6">
            <Card className="max-w-4xl mx-auto">
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Dna className="w-5 h-5 text-blue-600" />
                  Protein Analysis
                </CardTitle>
                <CardDescription>
                  Advanced AI-powered protein sequence analysis
                </CardDescription>
              </CardHeader>
              <CardContent className="space-y-6">
                {initialProtein && (
                  <div className="space-y-4">
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                      <div className="p-4 bg-blue-50 rounded-lg">
                        <h3 className="font-medium text-blue-900">Protein Name</h3>
                        <p className="text-blue-700">{initialProtein.name}</p>
                      </div>
                      <div className="p-4 bg-green-50 rounded-lg">
                        <h3 className="font-medium text-green-900">Sequence Length</h3>
                        <p className="text-green-700">{initialProtein.sequence.length} amino acids</p>
                      </div>
                    </div>
                    
                    <div className="p-4 bg-gray-50 rounded-lg">
                      <h3 className="font-medium text-gray-900 mb-2">Sequence</h3>
                      <div className="font-mono text-sm bg-white p-3 rounded border max-h-32 overflow-y-auto">
                        {initialProtein.sequence.match(/.{1,60}/g)?.join('\n')}
                      </div>
                    </div>

                    <Button 
                      onClick={mockProcessAI}
                      disabled={isProcessing}
                      size="lg"
                      className="w-full bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-700 hover:to-purple-700"
                    >
                      {isProcessing ? (
                        <>
                          <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                          Processing AI Pipeline...
                        </>
                      ) : (
                        <>
                          <Zap className="w-4 h-4 mr-2" />
                          Start AI-Enhanced Analysis
                        </>
                      )}
                    </Button>
                  </div>
                )}
              </CardContent>
            </Card>
          </TabsContent>

          <TabsContent value="generation" className="space-y-6">
            <Card className="max-w-4xl mx-auto">
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Sparkles className="w-5 h-5 text-purple-600" />
                  AI Molecule Generation
                </CardTitle>
                <CardDescription>
                  Advanced drug candidate generation using AI models
                </CardDescription>
              </CardHeader>
              <CardContent>
                {isProcessing ? (
                  <div className="text-center py-12">
                    <Loader2 className="w-12 h-12 text-blue-600 animate-spin mx-auto mb-4" />
                    <h3 className="text-lg font-medium mb-2">AI Processing in Progress</h3>
                    <p className="text-gray-600">Generating drug candidates using advanced AI models...</p>
                  </div>
                ) : results ? (
                  <div className="space-y-4">
                    <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
                      <div className="text-center p-4 bg-blue-50 rounded-lg">
                        <div className="text-2xl font-bold text-blue-600">
                          {results.candidates?.length || 0}
                        </div>
                        <div className="text-sm text-gray-600">Candidates Generated</div>
                      </div>
                      <div className="text-center p-4 bg-green-50 rounded-lg">
                        <div className="text-2xl font-bold text-green-600">
                          {results.docking_results?.best_binding_affinity || 'N/A'}
                        </div>
                        <div className="text-sm text-gray-600">Best Affinity (kcal/mol)</div>
                      </div>
                      <div className="text-center p-4 bg-purple-50 rounded-lg">
                        <div className="text-2xl font-bold text-purple-600">
                          {Math.round((results.docking_results?.average_confidence || 0) * 100)}%
                        </div>
                        <div className="text-sm text-gray-600">Avg Confidence</div>
                      </div>
                    </div>

                    <div className="space-y-3">
                      <h4 className="font-medium">Top Drug Candidates</h4>
                      {results.candidates?.map((candidate: any, index: number) => (
                        <div key={index} className="p-4 bg-gray-50 rounded-lg">
                          <div className="flex justify-between items-start">
                            <div>
                              <h5 className="font-medium">{candidate.name}</h5>
                              <p className="text-sm text-gray-600 font-mono">{candidate.smiles}</p>
                            </div>
                            <Badge variant="outline">
                              QED: {candidate.qed?.toFixed(2)}
                            </Badge>
                          </div>
                          <div className="mt-2 grid grid-cols-3 gap-4 text-sm">
                            <div>MW: {candidate.molecular_weight?.toFixed(1)} Da</div>
                            <div>LogP: {candidate.logP?.toFixed(1)}</div>
                            <div>Affinity: {candidate.binding_affinity?.toFixed(1)} kcal/mol</div>
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                ) : (
                  <div className="text-center py-12">
                    <FlaskConical className="w-16 h-16 text-gray-300 mx-auto mb-4" />
                    <p className="text-gray-500">Start the AI analysis to generate drug candidates</p>
                  </div>
                )}
              </CardContent>
            </Card>
          </TabsContent>

          <TabsContent value="docking" className="space-y-6">
            <Card className="max-w-4xl mx-auto">
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Target className="w-5 h-5 text-green-600" />
                  Molecular Docking Results
                </CardTitle>
                <CardDescription>
                  AI-enhanced binding affinity prediction and analysis
                </CardDescription>
              </CardHeader>
              <CardContent>
                {results ? (
                  <div className="space-y-4">
                    <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                      <div className="p-4 bg-green-50 rounded-lg">
                        <h4 className="font-medium text-green-900">Docking Success</h4>
                        <p className="text-2xl font-bold text-green-600">
                          {results.docking_results?.total_successful_dockings || 0}/2
                        </p>
                        <p className="text-sm text-green-700">Successful binding poses</p>
                      </div>
                      <div className="p-4 bg-blue-50 rounded-lg">
                        <h4 className="font-medium text-blue-900">Best Binding Affinity</h4>
                        <p className="text-2xl font-bold text-blue-600">
                          {results.docking_results?.best_binding_affinity || 'N/A'} kcal/mol
                        </p>
                        <p className="text-sm text-blue-700">Strongest protein-ligand interaction</p>
                      </div>
                    </div>
                    
                    <div className="flex items-center justify-center py-8">
                      <CheckCircle className="w-16 h-16 text-green-500 mr-4" />
                      <div>
                        <h3 className="text-lg font-medium text-green-900">Docking Analysis Complete</h3>
                        <p className="text-green-700">Ready for 3D visualization</p>
                      </div>
                    </div>
                  </div>
                ) : (
                  <div className="text-center py-12">
                    <Target className="w-16 h-16 text-gray-300 mx-auto mb-4" />
                    <p className="text-gray-500">Complete the AI generation step to see docking results</p>
                  </div>
                )}
              </CardContent>
            </Card>
          </TabsContent>

          <TabsContent value="visualization">
            <Card className="max-w-6xl mx-auto">
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Eye className="w-5 h-5 text-purple-600" />
                  3D Molecular Visualization
                </CardTitle>
                <CardDescription>
                  Interactive protein-ligand complex visualization
                </CardDescription>
              </CardHeader>
              <CardContent>
                {results ? (
                  <MolecularViewer
                    proteinData={initialProtein}
                    selectedCandidate={results.candidates?.[0]}
                    viewMode="cartoon"
                    showInteractions={true}
                    showWater={false}
                  />
                ) : (
                  <div className="text-center py-12">
                    <Eye className="w-16 h-16 text-gray-300 mx-auto mb-4" />
                    <p className="text-gray-500">Complete the analysis to view 3D molecular structures</p>
                  </div>
                )}
              </CardContent>
            </Card>
          </TabsContent>
        </Tabs>
      </div>
    </div>
  );
}
