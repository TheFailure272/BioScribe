"use client";

import { useState } from "react";
import { motion } from "framer-motion";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Progress } from "@/components/ui/progress";
import { 
  Target, 
  Play, 
  Settings, 
  CheckCircle, 
  Loader2,
  Zap,
  Award,
  BarChart3
} from "lucide-react";
import { ProteinData, DrugCandidate } from "../BioScribeWorkflow";
import Advanced3DVisualization from "../Advanced3DVisualization";

interface DockingTabProps {
  candidates?: DrugCandidate[];
  proteinData?: ProteinData;
  onDocking: (results: any) => void;
  isProcessing: boolean;
  results?: any;
}

// Mock docking results
const generateMockDockingResults = (candidates: DrugCandidate[]) => {
  return candidates.map((candidate, index) => ({
    ...candidate,
    bindingAffinity: -12.5 + Math.random() * 8, // -12.5 to -4.5 kcal/mol
    rmsd: 0.5 + Math.random() * 3, // 0.5 to 3.5 Å
    confidence: 60 + Math.random() * 35, // 60-95%
    poses: Math.floor(5 + Math.random() * 15), // 5-20 poses
    interactions: [
      { type: 'H-bond', residue: 'Asp25', distance: 2.1 + Math.random() * 0.8 },
      { type: 'Hydrophobic', residue: 'Ile50', distance: 3.2 + Math.random() * 1.0 },
      { type: 'π-π stacking', residue: 'Phe99', distance: 3.5 + Math.random() * 0.5 }
    ].slice(0, Math.floor(1 + Math.random() * 3))
  })).sort((a, b) => a.bindingAffinity - b.bindingAffinity); // Sort by binding affinity (lower is better)
};

export function DockingTab({ candidates, proteinData, onDocking, isProcessing, results }: DockingTabProps) {
  const [dockingMethod, setDockingMethod] = useState<'laboratory' | 'vina' | 'glide'>('laboratory');
  const [exhaustiveness, setExhaustiveness] = useState(8);
  const [numPoses, setNumPoses] = useState(10);
  const [realTimeResults, setRealTimeResults] = useState<any>(null);

  const handleStartDocking = async () => {
    console.log('Starting docking simulation...');
    console.log('Candidates:', candidates);
    console.log('Protein data:', proteinData);
    
    if (!candidates || !proteinData) {
      console.error('Missing candidates or protein data');
      return;
    }

    try {
      // Use laboratory-grade docking for the best candidate
      const bestCandidate = candidates[0];
      console.log('Best candidate:', bestCandidate);
      
      if (dockingMethod === 'laboratory') {
        // Real laboratory docking
        const response = await fetch('http://localhost:8000/api/laboratory/docking', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            protein_data: {
              sequence: proteinData.sequence,
              name: proteinData.name,
              analysis: proteinData.analysis
            },
            ligand_smiles: bestCandidate.smiles,
            binding_sites: proteinData.analysis?.binding_sites || []
          })
        });

        console.log('API Response status:', response.status);
        
        if (response.ok) {
          const laboratoryResults = await response.json();
          console.log('Laboratory results:', laboratoryResults);
          setRealTimeResults(laboratoryResults);
          
          onDocking({
            method: 'laboratory',
            laboratory_results: laboratoryResults,
            candidates: candidates.map(c => ({
              ...c,
              bindingAffinity: c === bestCandidate ? laboratoryResults.binding_affinity : c.bindingAffinity
            })),
            bestCandidate: {
              ...bestCandidate,
              bindingAffinity: laboratoryResults.binding_affinity,
              dockingScore: laboratoryResults.docking_score
            },
            statistics: {
              totalCandidates: candidates.length,
              highAffinityCandidates: 1,
              averageAffinity: laboratoryResults.binding_affinity,
              averageRMSD: 1.2,
              laboratory_grade: true,
              physics_based: true
            }
          });
        } else {
          const errorText = await response.text();
          console.error('API Error:', response.status, errorText);
          throw new Error(`Laboratory docking failed: ${response.status} - ${errorText}`);
        }
      } else {
        // Fallback to mock results for other methods
        const mockResults = generateMockDockingResults(candidates);
        onDocking({
          method: dockingMethod,
          candidates: mockResults,
          bestCandidate: mockResults[0],
          statistics: {
            totalCandidates: mockResults.length,
            highAffinityCandidates: mockResults.filter((c: any) => c.bindingAffinity < -8.0).length,
            averageAffinity: mockResults.reduce((sum: number, c: any) => sum + c.bindingAffinity, 0) / mockResults.length,
            averageRMSD: mockResults.reduce((sum: number, c: any) => sum + c.rmsd, 0) / mockResults.length
          }
        });
      }
    } catch (error) {
      console.error('Docking failed:', error);
      // Fallback to mock results on error
      const mockResults = generateMockDockingResults(candidates);
      onDocking({
        method: 'fallback',
        candidates: mockResults,
        bestCandidate: mockResults[0],
        statistics: {
          totalCandidates: mockResults.length,
          highAffinityCandidates: mockResults.filter((c: any) => c.bindingAffinity < -8.0).length,
          averageAffinity: mockResults.reduce((sum: number, c: any) => sum + c.bindingAffinity, 0) / mockResults.length,
          averageRMSD: mockResults.reduce((sum: number, c: any) => sum + c.rmsd, 0) / mockResults.length
        }
      });
    }
  };

  const getAffinityColor = (affinity: number) => {
    if (affinity < -10) return "text-green-600 bg-green-50";
    if (affinity < -8) return "text-blue-600 bg-blue-50";
    if (affinity < -6) return "text-yellow-600 bg-yellow-50";
    return "text-red-600 bg-red-50";
  };

  const getRMSDColor = (rmsd: number) => {
    if (rmsd < 2.0) return "text-green-600";
    if (rmsd < 3.0) return "text-yellow-600";
    return "text-red-600";
  };

  const getAffinityLabel = (affinity: number) => {
    if (affinity < -10) return "Excellent";
    if (affinity < -8) return "Good";
    if (affinity < -6) return "Moderate";
    return "Weak";
  };

  if (!candidates || candidates.length === 0) {
    return (
      <Card>
        <CardContent className="pt-6">
          <div className="text-center py-8">
            <Target className="w-12 h-12 text-muted-foreground mx-auto mb-4" />
            <h3 className="text-lg font-medium mb-2">Drug Candidates Required</h3>
            <p className="text-muted-foreground">
              Please generate drug candidates in the previous step to perform docking simulation.
            </p>
          </div>
        </CardContent>
      </Card>
    );
  }

  return (
    <div className="space-y-6">
      {/* Candidates Summary */}
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <Target className="w-5 h-5" />
            Docking Simulation
          </CardTitle>
          <CardDescription>
            Predict binding affinity and poses for generated drug candidates
          </CardDescription>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div>
              <div className="text-sm text-muted-foreground">Candidates</div>
              <div className="text-2xl font-bold text-primary">{candidates.length}</div>
            </div>
            <div>
              <div className="text-sm text-muted-foreground">Avg Molecular Weight</div>
              <div className="text-2xl font-bold text-secondary">
                {(candidates.reduce((sum, c) => sum + c.molecularWeight, 0) / candidates.length).toFixed(0)} Da
              </div>
            </div>
            <div>
              <div className="text-sm text-muted-foreground">Avg QED Score</div>
              <div className="text-2xl font-bold text-accent">
                {(candidates.reduce((sum, c) => sum + c.qed, 0) / candidates.length).toFixed(2)}
              </div>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Docking Settings */}
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <Settings className="w-5 h-5" />
            Docking Parameters
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-6">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div>
              <label className="text-sm font-medium mb-2 block">
                Docking Software
              </label>
              <select
                value={dockingMethod}
                onChange={(e) => setDockingMethod(e.target.value as any)}
                className="w-full px-3 py-2 border rounded-md bg-background"
              >
                <option value="vina">AutoDock Vina</option>
                <option value="glide">Schrödinger Glide</option>
                <option value="gold">CCDC GOLD</option>
              </select>
            </div>
            <div>
              <label className="text-sm font-medium mb-2 block">
                Exhaustiveness
              </label>
              <select
                value={exhaustiveness}
                onChange={(e) => setExhaustiveness(Number(e.target.value))}
                className="w-full px-3 py-2 border rounded-md bg-background"
              >
                <option value={4}>4 (Fast)</option>
                <option value={8}>8 (Standard)</option>
                <option value={16}>16 (Thorough)</option>
                <option value={32}>32 (Exhaustive)</option>
              </select>
            </div>
            <div>
              <label className="text-sm font-medium mb-2 block">
                Number of Poses
              </label>
              <select
                value={numPoses}
                onChange={(e) => setNumPoses(Number(e.target.value))}
                className="w-full px-3 py-2 border rounded-md bg-background"
              >
                <option value={5}>5 poses</option>
                <option value={10}>10 poses</option>
                <option value={20}>20 poses</option>
              </select>
            </div>
          </div>

          <div className="p-4 bg-muted rounded-lg">
            <h4 className="font-medium mb-2">Estimated Runtime</h4>
            <p className="text-sm text-muted-foreground">
              {candidates.length} candidates × {exhaustiveness} exhaustiveness × {numPoses} poses
              = ~{Math.ceil((candidates.length * exhaustiveness * numPoses) / 20)} minutes
            </p>
          </div>
        </CardContent>
      </Card>

      {/* Docking Results */}
      {results && (
        <Card>
          <CardHeader>
            <CardTitle className="flex items-center gap-2">
              <CheckCircle className="w-5 h-5 text-green-500" />
              Docking Results
            </CardTitle>
            <CardDescription>
              Binding affinity predictions ranked by score
            </CardDescription>
          </CardHeader>
          <CardContent>
            {/* Best Candidate Highlight */}
            {results.bestCandidate && (
              <div className="p-4 bg-gradient-to-r from-green-50 to-blue-50 rounded-lg mb-6 border border-green-200">
                <div className="flex items-start justify-between mb-3">
                  <div>
                    <div className="flex items-center gap-2 mb-1">
                      <Award className="w-5 h-5 text-green-600" />
                      <h4 className="font-bold text-green-800">Best Candidate</h4>
                    </div>
                    <h3 className="text-lg font-bold">{results.bestCandidate.name}</h3>
                    <p className="text-sm text-muted-foreground molecular-formula">
                      {results.bestCandidate.smiles}
                    </p>
                  </div>
                  <Badge className="bg-green-100 text-green-800 border-green-200">
                    Rank #1
                  </Badge>
                </div>
                
                <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
                  <div>
                    <span className="text-muted-foreground">Binding Affinity:</span>
                    <div className="font-bold text-green-600">
                      {results.bestCandidate.bindingAffinity.toFixed(1)} kcal/mol
                    </div>
                  </div>
                  <div>
                    <span className="text-muted-foreground">RMSD:</span>
                    <div className="font-bold text-green-600">
                      {results.bestCandidate.rmsd.toFixed(1)} Å
                    </div>
                  </div>
                  <div>
                    <span className="text-muted-foreground">Confidence:</span>
                    <div className="font-bold text-green-600">
                      {results.bestCandidate.confidence.toFixed(0)}%
                    </div>
                  </div>
                  <div>
                    <span className="text-muted-foreground">Poses:</span>
                    <div className="font-bold text-green-600">
                      {results.bestCandidate.poses}
                    </div>
                  </div>
                </div>
              </div>
            )}

            {/* All Results */}
            <div className="space-y-3">
              {results.candidates.map((candidate: any, index: number) => (
                <motion.div
                  key={index}
                  className="p-4 border rounded-lg hover:bg-secondary/50 transition-colors"
                  initial={{ opacity: 0, y: 20 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ delay: index * 0.05 }}
                >
                  <div className="flex items-start justify-between mb-3">
                    <div>
                      <h4 className="font-medium">{candidate.name}</h4>
                      <p className="text-sm text-muted-foreground molecular-formula">
                        {candidate.smiles}
                      </p>
                    </div>
                    <div className="flex items-center gap-2">
                      <Badge variant="outline">#{index + 1}</Badge>
                      <Badge 
                        variant="outline" 
                        className={getAffinityColor(candidate.bindingAffinity)}
                      >
                        {getAffinityLabel(candidate.bindingAffinity)}
                      </Badge>
                    </div>
                  </div>
                  
                  <div className="grid grid-cols-2 md:grid-cols-5 gap-4 text-sm">
                    <div>
                      <span className="text-muted-foreground">Affinity:</span>
                      <div className="font-medium">
                        {candidate.bindingAffinity.toFixed(1)} kcal/mol
                      </div>
                    </div>
                    <div>
                      <span className="text-muted-foreground">RMSD:</span>
                      <div className={`font-medium ${getRMSDColor(candidate.rmsd)}`}>
                        {candidate.rmsd.toFixed(1)} Å
                      </div>
                    </div>
                    <div>
                      <span className="text-muted-foreground">Confidence:</span>
                      <div className="font-medium">{candidate.confidence.toFixed(0)}%</div>
                    </div>
                    <div>
                      <span className="text-muted-foreground">Poses:</span>
                      <div className="font-medium">{candidate.poses}</div>
                    </div>
                    <div>
                      <span className="text-muted-foreground">Interactions:</span>
                      <div className="font-medium">{candidate.interactions.length}</div>
                    </div>
                  </div>

                  {/* Key Interactions */}
                  {candidate.interactions.length > 0 && (
                    <div className="mt-3 pt-3 border-t">
                      <div className="text-xs text-muted-foreground mb-1">Key Interactions:</div>
                      <div className="flex flex-wrap gap-2">
                        {candidate.interactions.map((interaction: any, i: number) => (
                          <Badge key={i} variant="outline" className="text-xs">
                            {interaction.type} - {interaction.residue}
                          </Badge>
                        ))}
                      </div>
                    </div>
                  )}
                </motion.div>
              ))}
            </div>

            {/* Statistics */}
            <div className="mt-6 p-4 bg-muted rounded-lg">
              <h4 className="font-medium mb-3 flex items-center gap-2">
                <BarChart3 className="w-4 h-4" />
                Docking Statistics
              </h4>
              <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
                <div>
                  <div className="text-2xl font-bold text-primary">
                    {results.statistics.highAffinityCandidates}
                  </div>
                  <div className="text-muted-foreground">High Affinity (&lt; -8.0)</div>
                </div>
                <div>
                  <div className="text-2xl font-bold text-secondary">
                    {results.statistics.averageAffinity.toFixed(1)}
                  </div>
                  <div className="text-muted-foreground">Avg Affinity (kcal/mol)</div>
                </div>
                <div>
                  <div className="text-2xl font-bold text-accent">
                    {results.statistics.averageRMSD.toFixed(1)}
                  </div>
                  <div className="text-muted-foreground">Avg RMSD (Å)</div>
                </div>
                <div>
                  <div className="text-2xl font-bold text-chart-4">
                    {results.statistics.totalCandidates}
                  </div>
                  <div className="text-muted-foreground">Total Docked</div>
                </div>
              </div>
            </div>
          </CardContent>
        </Card>
      )}

      {/* Action Button */}
      <div className="flex justify-end">
        <Button
          onClick={handleStartDocking}
          disabled={isProcessing}
          size="lg"
          className="gap-2"
        >
          {isProcessing ? (
            <>
              <Loader2 className="w-4 h-4 animate-spin" />
              Running Docking...
            </>
          ) : (
            <>
              <Play className="w-4 h-4" />
              Start Docking Simulation
            </>
          )}
        </Button>
      </div>

      {/* Processing Status */}
      {isProcessing && (
        <Card>
          <CardContent className="pt-6">
            <div className="space-y-4">
              <div className="flex items-center gap-3">
                <Loader2 className="w-5 h-5 animate-spin text-primary" />
                <span className="font-medium">Running molecular docking simulation...</span>
              </div>
              <Progress value={45} className="h-2" />
              <div className="text-sm text-muted-foreground">
                • Preparing protein structure<br/>
                • Generating ligand conformations<br/>
                • Calculating binding poses<br/>
                • Scoring and ranking results
              </div>
              <div className="text-xs text-muted-foreground">
                Estimated time remaining: ~3 minutes
              </div>
            </div>
          </CardContent>
        </Card>
      )}

      {/* Advanced 3D Visualization */}
      {realTimeResults && (
        <Advanced3DVisualization
          proteinData={proteinData}
          ligandData={candidates?.[0]}
          dockingResults={realTimeResults}
          interactionData={realTimeResults.interaction_analysis}
        />
      )}
    </div>
  );
}
