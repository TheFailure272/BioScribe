'use client';

import React, { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import { ArrowLeft, Download, Share2, Sparkles } from 'lucide-react';
import Link from 'next/link';
import { ExecutiveResults } from '@/components/ExecutiveResults';

export default function ResultsPage() {
    const [isLoading, setIsLoading] = useState(true);

    useEffect(() => {
        setTimeout(() => setIsLoading(false), 2000);
    }, []);

    const mockResults = {
        results: {
            overall_executive_summary: {
                title: "Drug Discovery Pipeline - Complete Analysis",
                subtitle: "Target: Protein Kinase (1ATP)",
                completion_time: "2 minutes",
                key_findings: [
                    "20 high-affinity drug candidates identified",
                    "Best binding affinity: -9.2 kcal/mol",
                    "3 candidates pass all ADMET criteria",
                    "Synthesis routes validated"
                ],
                key_achievements: [
                    "Successfully analyzed protein structure",
                    "Generated 20 novel drug candidates",
                    "Validated binding affinities",
                    "Confirmed synthesis feasibility"
                ],
                recommendations: [
                    "Proceed with top 3 candidates for in vitro testing",
                    "Validate synthesis routes experimentally",
                    "Consider structural modifications for improved selectivity",
                    "Initiate DMPK studies for lead optimization"
                ],
                next_steps: [
                    "1. Conduct in vitro binding assays for top candidates",
                    "2. Perform toxicity screening on lead compounds",
                    "3. Optimize synthesis routes for scale-up",
                    "4. Begin pharmacokinetic studies"
                ],
                pipeline_statistics: {
                    total_steps_completed: 6,
                    total_time_seconds: 120,
                    success_rate: 100
                }
            },
            protein_analysis_summary: {
                key_findings: [
                    "3 binding sites identified",
                    "Drugability score: 0.85",
                    "Active site volume: 450 Å³"
                ]
            },
            drug_generation_summary: {
                key_findings: [
                    "20 candidates generated",
                    "Average MW: 342.5 Da",
                    "85% drug-like"
                ]
            },
            docking_summary: {
                key_findings: [
                    "Best affinity: -9.2 kcal/mol",
                    "12 key interactions",
                    "100% success rate"
                ]
            },
            blockchain_summary: {
                key_findings: [
                    "Data stored on Ethereum",
                    "IPFS hash: QmXyz123",
                    "Immutable record created"
                ]
            },
            fair_summary: {
                key_findings: [
                    "100% FAIR compliant",
                    "All data findable",
                    "Fully interoperable"
                ]
            }
        }
    };

    if (isLoading) {
        return (
            <div className="min-h-screen bg-gradient-to-br from-violet-100 via-purple-50 to-fuchsia-100 flex items-center justify-center">
                <div className="text-center">
                    <motion.div
                        animate={{ rotate: 360 }}
                        transition={{ duration: 2, repeat: Infinity, ease: "linear" }}
                        className="w-20 h-20 rounded-full bg-gradient-to-br from-violet-600 to-fuchsia-600 flex items-center justify-center mx-auto mb-6 shadow-2xl shadow-violet-500/50"
                    >
                        < Sparkles className="w-10 h-10 text-white" />
                    </motion.div>
                    <h2 className="text-3xl font-light text-slate-900 mb-2">Analyzing...</h2>
                    <p className="text-lg font-light text-slate-600">Generating insights</p>
                </div>
            </div>
        );
    }

    return (
        <div className="min-h-screen bg-gradient-to-br from-violet-100 via-purple-50 to-fuchsia-100">
            <div className="bg-white/70 backdrop-blur-xl border-b border-slate-100/50 sticky top-20 z-40">
                <div className="max-w-7xl mx-auto px-8 py-6">
                    <div className="flex items-center justify-between">
                        <div className="flex items-center gap-6">
                            <Link href="/dashboard">
                                <button className="p-2 rounded-xl hover:bg-slate-100 transition-colors">
                                    <ArrowLeft className="w-5 h-5 text-slate-600" />
                                </button>
                            </Link>
                            <div>
                                <h1 className="text-3xl font-light text-slate-900">Pipeline Results</h1>
                                <p className="text-sm font-light text-slate-600">Target: 1ATP • Completed 2 min ago</p>
                            </div>
                        </div>
                        <div className="flex items-center gap-3">
                            <button className="px-5 py-2.5 rounded-xl border-2 border-slate-200 text-slate-600 font-light hover:bg-white/50 transition-all flex items-center gap-2">
                                <Share2 className="w-4 h-4" />
                                Share
                            </button>
                            <button className="px-5 py-2.5 rounded-xl bg-gradient-to-r from-violet-600 to-fuchsia-600 text-white font-light shadow-lg shadow-violet-500/30 hover:shadow-xl transition-all flex items-center gap-2">
                                <Download className="w-4 h-4" />
                                Export
                            </button>
                        </div>
                    </div>
                </div>
            </div>
            <div className="max-w-7xl mx-auto px-8 py-12">
                <ExecutiveResults results={mockResults} />
            </div>
        </div>
    );
}
