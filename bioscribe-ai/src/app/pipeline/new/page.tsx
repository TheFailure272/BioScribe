'use client';

import React, { useState } from 'react';
import { useRouter } from 'next/navigation';
import {
    Upload,
    Settings,
    Play,
    CheckCircle2,
    ArrowRight,
    ArrowLeft,
    Sparkles
} from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';

export default function PipelineWizardPage() {
    const router = useRouter();
    const [currentStep, setCurrentStep] = useState(1);
    const [isRunning, setIsRunning] = useState(false);
    const [formData, setFormData] = useState({
        proteinFile: null as File | null,
        pdbId: '',
        generateMolecules: true,
        admetAnalysis: true,
        retrosynthesis: true
    });

    const steps = [
        { number: 1, title: 'Upload', subtitle: 'Target protein' },
        { number: 2, title: 'Configure', subtitle: 'Analysis options' },
        { number: 3, title: 'Run', subtitle: 'Start discovery' }
    ];

    const handleNext = () => {
        if (currentStep < 3) {
            setCurrentStep(currentStep + 1);
        } else {
            setIsRunning(true);
            setTimeout(() => {
                router.push('/results/1');
            }, 3500);
        }
    };

    const handleBack = () => {
        if (currentStep > 1) {
            setCurrentStep(currentStep - 1);
        }
    };

    return (
        <div className="min-h-screen bg-gradient-to-br from-violet-100 via-purple-50 to-fuchsia-100 py-24">
            <div className="max-w-4xl mx-auto px-8">
                {/* Minimal Progress Dots */}
                <div className="mb-16">
                    <div className="flex items-center justify-center gap-3">
                        {steps.map((step, idx) => (
                            <React.Fragment key={step.number}>
                                <motion.div className="flex flex-col items-center gap-3">
                                    <motion.div
                                        animate={{
                                            scale: currentStep === step.number ? 1 : 0.7,
                                            opacity: currentStep >= step.number ? 1 : 0.3
                                        }}
                                        transition={{ duration: 0.5, ease: "easeOut" }}
                                        className={`w-3 h-3 rounded-full ${currentStep >= step.number
                                                ? 'bg-gradient-to-r from-violet-600 to-fuchsia-600'
                                                : 'bg-slate-300'
                                            }`}
                                    />
                                    <div className="text-center">
                                        <div className={`text-sm font-light ${currentStep === step.number ? 'text-violet-600' : 'text-slate-400'
                                            }`}>
                                            {step.title}
                                        </div>
                                        <div className="text-xs font-light text-slate-400">{step.subtitle}</div>
                                    </div>
                                </motion.div>
                                {idx < steps.length - 1 && (
                                    <div className="w-16 h-px bg-slate-200 mb-8" />
                                )}
                            </React.Fragment>
                        ))}
                    </div>
                </div>

                {/* Step Content Card */}
                <AnimatePresence mode="wait">
                    <motion.div
                        key={currentStep}
                        initial={{ opacity: 0, scale: 0.98 }}
                        animate={{ opacity: 1, scale: 1 }}
                        exit={{ opacity: 0, scale: 0.98 }}
                        transition={{ duration: 0.5, ease: "easeOut" }}
                        className="bg-white/60 backdrop-blur-sm rounded-3xl shadow-2xl p-12 mb-12 border border-slate-100/50"
                    >
                        {currentStep === 1 && <Step1 formData={formData} setFormData={setFormData} />}
                        {currentStep === 2 && <Step2 formData={formData} setFormData={setFormData} />}
                        {currentStep === 3 && <Step3 isRunning={isRunning} />}
                    </motion.div>
                </AnimatePresence>

                {/* Navigation */}
                <div className="flex justify-between items-center">
                    <motion.button
                        onClick={handleBack}
                        disabled={currentStep === 1}
                        whileHover={{ scale: 1.02 }}
                        whileTap={{ scale: 0.98 }}
                        className="px-8 py-4 rounded-2xl border-2 border-slate-200 text-slate-600 font-light hover:bg-white/50 hover:border-slate-300 disabled:opacity-30 disabled:cursor-not-allowed transition-all duration-400 flex items-center gap-2"
                    >
                        <ArrowLeft className="w-5 h-5" />
                        Back
                    </motion.button>

                    <motion.button
                        onClick={handleNext}
                        disabled={isRunning}
                        whileHover={{ scale: 1.02 }}
                        whileTap={{ scale: 0.98 }}
                        className="px-10 py-4 rounded-2xl bg-gradient-to-r from-violet-600 to-fuchsia-600 text-white font-light shadow-2xl shadow-violet-500/50 hover:shadow-violet-500/70 disabled:opacity-50 disabled:cursor-not-allowed transition-all duration-400 flex items-center gap-3"
                    >
                        {currentStep === 3 ? (
                            <>
                                {isRunning ? 'Running...' : (
                                    <>
                                        <Play className="w-5 h-5" />
                                        Run Pipeline
                                    </>
                                )}
                            </>
                        ) : (
                            <>
                                Continue
                                <ArrowRight className="w-5 h-5" />
                            </>
                        )}
                    </motion.button>
                </div>
            </div>
        </div>
    );
}

function Step1({ formData, setFormData }: any) {
    const examples = [
        { id: '1HIV', name: 'HIV-1 Protease', description: 'HIV drug target' },
        { id: '6LU7', name: 'SARS-CoV-2 Mpro', description: 'COVID-19 target' },
        { id: '1M17', name: 'EGFR Kinase', description: 'Cancer target' },
        { id: '1ATP', name: 'Protein Kinase', description: 'ATP binding' }
    ];

    return (
        <div>
            <h2 className="text-5xl font-light text-slate-900 mb-3">Upload Target</h2>
            <p className="text-xl font-light text-slate-600 mb-12">Begin with your protein structure</p>

            {/* File Upload Area */}
            <div className="mb-12">
                <label className="block text-sm font-light text-slate-500 mb-4 uppercase tracking-wider">PDB File</label>
                <div className="relative group">
                    <div className="absolute inset-0 bg-gradient-to-r from-violet-600 to-fuchsia-600 rounded-2xl blur-xl opacity-20 group-hover:opacity-30 transition-opacity duration-500" />
                    <div className="relative border-2 border-dashed border-slate-200 group-hover:border-violet-300 rounded-2xl p-16 text-center transition-all duration-500 bg-white/40 backdrop-blur-sm cursor-pointer">
                        <Upload className="w-12 h-12 text-slate-300 mx-auto mb-4" />
                        <p className="text-slate-600 font-light text-lg">Drop your file here</p>
                        <p className="text-sm font-light text-slate-400 mt-2">or click to browse</p>
                        <input type="file" accept=".pdb" className="hidden" />
                    </div>
                </div>
            </div>

            <div className="flex items-center gap-6 my-12">
                <div className="flex-1 h-px bg-slate-200" />
                <span className="text-sm font-light text-slate-400">OR</span>
                <div className="flex-1 h-px bg-slate-200" />
            </div>

            {/* PDB ID Input */}
            <div className="mb-12">
                <label className="block text-sm font-light text-slate-500 mb-4 uppercase tracking-wider">PDB ID</label>
                <input
                    type="text"
                    placeholder="Enter PDB ID (e.g., 1ATP)"
                    value={formData.pdbId}
                    onChange={(e) => setFormData({ ...formData, pdbId: e.target.value })}
                    className="w-full px-6 py-4 rounded-xl border border-slate-200 focus:border-violet-300 focus:ring-4 focus:ring-violet-100 outline-none transition-all duration-400 text-lg font-light bg-white/40 backdrop-blur-sm"
                />
                <p className="text-sm font-light text-slate-400 mt-3">Fetched from Protein Data Bank</p>
            </div>

            {/* Example Proteins */}
            <div>
                <label className="block text-sm font-light text-slate-500 mb-4 uppercase tracking-wider">Or Try Examples</label>
                <div className="grid grid-cols-2 gap-3">
                    {examples.map((example) => (
                        <button
                            key={example.id}
                            onClick={() => setFormData({ ...formData, pdbId: example.id })}
                            className="text-left p-4 rounded-xl border-2 border-slate-200 hover:border-violet-300 hover:bg-violet-50/30 transition-all duration-300"
                        >
                            <div className="font-normal text-slate-900">{example.name}</div>
                            <div className="text-sm font-light text-slate-500">{example.description}</div>
                            <div className="text-xs font-mono text-violet-600 mt-1">{example.id}</div>
                        </button>
                    ))}
                </div>
            </div>
        </div>
    );
}

function Step2({ formData, setFormData }: any) {
    const options = [
        { key: 'generateMolecules', label: 'De Novo Generation', description: 'AI-designed molecules' },
        { key: 'admetAnalysis', label: 'ADMET Prediction', description: 'Drug-likeness analysis' },
        { key: 'retrosynthesis', label: 'Retrosynthesis', description: 'Synthesis planning' }
    ];

    return (
        <div>
            <h2 className="text-5xl font-light text-slate-900 mb-3">Configure</h2>
            <p className="text-xl font-light text-slate-600 mb-12">Select your analyses</p>

            <div className="space-y-5">
                {options.map((option) => (
                    <label
                        key={option.key}
                        className="group flex items-center gap-5 p-6 rounded-2xl border-2 border-slate-200 hover:border-violet-300 hover:bg-violet-50/30 cursor-pointer transition-all duration-500"
                    >
                        <div className={`w-6 h-6 rounded-full border-2 flex items-center justify-center transition-all duration-400 ${formData[option.key]
                                ? 'border-violet-600 bg-gradient-to-r from-violet-600 to-fuchsia-600'
                                : 'border-slate-300 group-hover:border-violet-300'
                            }`}>
                            {formData[option.key] && <CheckCircle2 className="w-4 h-4 text-white" />}
                        </div>
                        <input
                            type="checkbox"
                            checked={formData[option.key]}
                            onChange={(e) => setFormData({ ...formData, [option.key]: e.target.checked })}
                            className="hidden"
                        />
                        <div className="flex-1">
                            <div className="font-normal text-lg text-slate-900">{option.label}</div>
                            <div className="text-sm font-light text-slate-500">{option.description}</div>
                        </div>
                    </label>
                ))}
            </div>
        </div>
    );
}

function Step3({ isRunning }: { isRunning: boolean }) {
    return (
        <div className="text-center py-16">
            {!isRunning ? (
                <>
                    <motion.div
                        initial={{ scale: 0 }}
                        animate={{ scale: 1 }}
                        transition={{ duration: 0.6, ease: "easeOut" }}
                        className="w-24 h-24 rounded-full bg-gradient-to-br from-emerald-500 to-green-500 flex items-center justify-center mx-auto mb-8 shadow-2xl shadow-emerald-500/50"
                    >
                        <CheckCircle2 className="w-12 h-12 text-white" />
                    </motion.div>
                    <h2 className="text-5xl font-light text-slate-900 mb-3">Ready</h2>
                    <p className="text-xl font-light text-slate-600 mb-12">Everything is configured</p>
                    <div className="bg-gradient-to-r from-violet-50 to-fuchsia-50 rounded-2xl p-8 max-w-md mx-auto border border-violet-100">
                        <div className="space-y-4 font-light">
                            <div className="flex justify-between text-slate-600">
                                <span>Target</span>
                                <span className="font-normal text-slate-900">PDB: 1ATP</span>
                            </div>
                            <div className="flex justify-between text-slate-600">
                                <span>Analyses</span>
                                <span className="font-normal text-slate-900">3 selected</span>
                            </div>
                            <div className="flex justify-between text-slate-600">
                                <span>Est. Time</span>
                                <span className="font-normal text-slate-900">~5 minutes</span>
                            </div>
                        </div>
                    </div>
                </>
            ) : (
                <>
                    <motion.div
                        animate={{ rotate: 360 }}
                        transition={{ duration: 3, repeat: Infinity, ease: "linear" }}
                        className="w-24 h-24 rounded-full bg-gradient-to-br from-violet-600 to-fuchsia-600 flex items-center justify-center mx-auto mb-8 shadow-2xl shadow-violet-500/50"
                    >
                        <Sparkles className="w-12 h-12 text-white" />
                    </motion.div>
                    <h2 className="text-5xl font-light text-slate-900 mb-3">Running</h2>
                    <p className="text-xl font-light text-slate-600 mb-12">Analyzing your protein</p>
                    <div className="max-w-md mx-auto">
                        <div className="h-1 bg-slate-200 rounded-full overflow-hidden">
                            <motion.div
                                initial={{ width: "0%" }}
                                animate={{ width: "100%" }}
                                transition={{ duration: 3.5, ease: "easeInOut" }}
                                className="h-full bg-gradient-to-r from-violet-600 to-fuchsia-600"
                            />
                        </div>
                    </div>
                </>
            )}
        </div>
    );
}
