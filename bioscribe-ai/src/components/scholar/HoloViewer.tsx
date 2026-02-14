'use client';

import React, { useState, useRef, useCallback, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import Webcam from 'react-webcam';
import { useGamificationStore } from '@/stores/GamificationStore';

// ─── Sample Molecules ────────────────────────────────────────────────
const MOLECULES = [
    {
        id: 'aspirin',
        name: 'Aspirin (Acetylsalicylic Acid)',
        formula: 'C₉H₈O₄',
        atoms: [
            { element: 'C', x: 0, y: 0, z: 0, color: '#333333' },
            { element: 'C', x: 1.2, y: 0.5, z: 0, color: '#333333' },
            { element: 'C', x: 2.4, y: 0, z: 0, color: '#333333' },
            { element: 'C', x: 2.4, y: -1.5, z: 0, color: '#333333' },
            { element: 'C', x: 1.2, y: -2, z: 0, color: '#333333' },
            { element: 'C', x: 0, y: -1.5, z: 0, color: '#333333' },
            { element: 'O', x: -1.2, y: 0.5, z: 0, color: '#dc2626' },
            { element: 'O', x: 3.6, y: 0.5, z: 0, color: '#dc2626' },
            { element: 'C', x: -1.2, y: 2, z: 0, color: '#333333' },
            { element: 'O', x: -2.4, y: 2.5, z: 0, color: '#dc2626' },
            { element: 'C', x: 4.8, y: 0, z: 0, color: '#333333' },
            { element: 'O', x: 4.8, y: -1.5, z: 0, color: '#dc2626' },
        ],
        bonds: [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6], [2, 7], [6, 8], [8, 9], [7, 10], [10, 11],
        ],
        description: 'Non-selective COX inhibitor — antiplatelet, anti-inflammatory, antipyretic',
        teachingPoint: '💡 Aspirin irreversibly inhibits COX-1, preventing thromboxane A2 synthesis in platelets for their entire 7-10 day lifespan.',
    },
    {
        id: 'caffeine',
        name: 'Caffeine',
        formula: 'C₈H₁₀N₄O₂',
        atoms: [
            { element: 'C', x: 0, y: 0, z: 0, color: '#333333' },
            { element: 'N', x: 1.2, y: 0.5, z: 0, color: '#2563eb' },
            { element: 'C', x: 2.4, y: 0, z: 0, color: '#333333' },
            { element: 'C', x: 2.4, y: -1.5, z: 0, color: '#333333' },
            { element: 'N', x: 1.2, y: -2, z: 0, color: '#2563eb' },
            { element: 'C', x: 0, y: -1.5, z: 0, color: '#333333' },
            { element: 'N', x: 3.6, y: 0.5, z: 0, color: '#2563eb' },
            { element: 'C', x: 3.6, y: -2, z: 0, color: '#333333' },
            { element: 'N', x: 4.5, y: -1, z: 0, color: '#2563eb' },
            { element: 'O', x: -1.2, y: 0.5, z: 0, color: '#dc2626' },
            { element: 'O', x: 0, y: -3, z: 0, color: '#dc2626' },
        ],
        bonds: [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [2, 6], [6, 8], [8, 7], [7, 3], [0, 9], [5, 10],
        ],
        description: 'Adenosine receptor antagonist — blocks drowsiness signals',
        teachingPoint: '💡 Caffeine blocks A1 and A2A adenosine receptors, preventing the sleep-promoting effects of adenosine accumulation.',
    },
    {
        id: 'penicillin',
        name: 'Penicillin G',
        formula: 'C₁₆H₁₈N₂O₄S',
        atoms: [
            { element: 'C', x: 0, y: 0, z: 0, color: '#333333' },
            { element: 'C', x: 1.2, y: 0.5, z: 0, color: '#333333' },
            { element: 'N', x: 2.0, y: -0.3, z: 0, color: '#2563eb' },
            { element: 'C', x: 1.5, y: -1.5, z: 0, color: '#333333' },
            { element: 'S', x: 0.3, y: -1, z: 0, color: '#eab308' },
            { element: 'C', x: 3.2, y: 0, z: 0, color: '#333333' },
            { element: 'O', x: 3.2, y: 1.3, z: 0, color: '#dc2626' },
            { element: 'N', x: 4.2, y: -0.5, z: 0, color: '#2563eb' },
            { element: 'O', x: -1, y: 0.8, z: 0, color: '#dc2626' },
            { element: 'O', x: -1, y: -0.8, z: 0, color: '#dc2626' },
        ],
        bonds: [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 0], [2, 5], [5, 6], [5, 7], [0, 8], [0, 9],
        ],
        description: 'β-lactam antibiotic — inhibits bacterial cell wall synthesis',
        teachingPoint: '💡 The β-lactam ring (the 4-membered ring with N) irreversibly binds transpeptidase (PBP), preventing crosslinking of peptidoglycan in bacterial cell walls.',
    },
];

export default function HoloViewer() {
    const [selectedMolecule, setSelectedMolecule] = useState(MOLECULES[0]);
    const [showWebcam, setShowWebcam] = useState(false);
    const [rotation, setRotation] = useState({ x: 0, y: 0 });
    const [scale, setScale] = useState(1);
    const [isDragging, setIsDragging] = useState(false);
    const [showInfo, setShowInfo] = useState(true);
    const lastPos = useRef({ x: 0, y: 0 });
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const animRef = useRef<number>(0);

    const { awardXP } = useGamificationStore();

    // ─── Canvas Rendering ──────────────────────────────────────────
    const renderMolecule = useCallback(() => {
        const canvas = canvasRef.current;
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        if (!ctx) return;

        const w = canvas.width;
        const h = canvas.height;
        ctx.clearRect(0, 0, w, h);

        const cx = w / 2;
        const cy = h / 2;
        const s = 40 * scale;

        const cosX = Math.cos(rotation.x);
        const sinX = Math.sin(rotation.x);
        const cosY = Math.cos(rotation.y);
        const sinY = Math.sin(rotation.y);

        // Project 3D → 2D
        const projected = selectedMolecule.atoms.map((atom) => {
            let x = atom.x, y = atom.y, z = atom.z;
            const y1 = y * cosX - z * sinX;
            const z1 = y * sinX + z * cosX;
            const x1 = x * cosY + z1 * sinY;
            return {
                sx: cx + x1 * s,
                sy: cy + y1 * s,
                z: z1,
                element: atom.element,
                color: atom.color,
            };
        });

        // Draw bonds
        ctx.lineWidth = 2;
        for (const [a, b] of selectedMolecule.bonds) {
            const pa = projected[a], pb = projected[b];
            if (!pa || !pb) continue;
            ctx.beginPath();
            ctx.moveTo(pa.sx, pa.sy);
            ctx.lineTo(pb.sx, pb.sy);
            ctx.strokeStyle = showWebcam ? 'rgba(0, 240, 255, 0.6)' : 'rgba(255, 255, 255, 0.3)';
            ctx.stroke();
        }

        // Draw atoms
        for (const p of projected) {
            const radius = p.element === 'H' ? 8 : p.element === 'C' ? 12 : 11;
            ctx.beginPath();
            ctx.arc(p.sx, p.sy, radius * scale, 0, Math.PI * 2);

            if (showWebcam) {
                ctx.fillStyle = p.color;
                ctx.globalAlpha = 0.85;
                ctx.fill();
                ctx.globalAlpha = 1;
                ctx.strokeStyle = '#00f0ff';
                ctx.lineWidth = 1.5;
                ctx.stroke();
                // Glow
                ctx.shadowColor = p.color;
                ctx.shadowBlur = 15;
                ctx.beginPath();
                ctx.arc(p.sx, p.sy, radius * scale, 0, Math.PI * 2);
                ctx.fill();
                ctx.shadowBlur = 0;
            } else {
                ctx.fillStyle = p.color;
                ctx.fill();
                ctx.strokeStyle = 'rgba(255,255,255,0.2)';
                ctx.lineWidth = 1;
                ctx.stroke();
            }

            // Element label
            ctx.fillStyle = showWebcam ? '#00f0ff' : '#fff';
            ctx.font = `${10 * scale}px JetBrains Mono, monospace`;
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText(p.element, p.sx, p.sy);
        }

        // Auto-rotate when not dragging
        if (!isDragging) {
            setRotation((prev) => ({ ...prev, y: prev.y + 0.005 }));
        }

        animRef.current = requestAnimationFrame(renderMolecule);
    }, [selectedMolecule, rotation, scale, showWebcam, isDragging]);

    useEffect(() => {
        animRef.current = requestAnimationFrame(renderMolecule);
        return () => cancelAnimationFrame(animRef.current);
    }, [renderMolecule]);

    // ─── Resize Canvas ─────────────────────────────────────────────
    useEffect(() => {
        const resize = () => {
            const canvas = canvasRef.current;
            if (!canvas) return;
            const parent = canvas.parentElement;
            if (!parent) return;
            canvas.width = parent.clientWidth;
            canvas.height = parent.clientHeight;
        };
        resize();
        window.addEventListener('resize', resize);
        return () => window.removeEventListener('resize', resize);
    }, []);

    // ─── Mouse Interaction ─────────────────────────────────────────
    const handleMouseDown = (e: React.MouseEvent) => {
        setIsDragging(true);
        lastPos.current = { x: e.clientX, y: e.clientY };
    };

    const handleMouseMove = (e: React.MouseEvent) => {
        if (!isDragging) return;
        const dx = e.clientX - lastPos.current.x;
        const dy = e.clientY - lastPos.current.y;
        setRotation((prev) => ({
            x: prev.x + dy * 0.01,
            y: prev.y + dx * 0.01,
        }));
        lastPos.current = { x: e.clientX, y: e.clientY };
    };

    const handleWheel = (e: React.WheelEvent) => {
        e.preventDefault();
        setScale((prev) => Math.max(0.3, Math.min(3, prev - e.deltaY * 0.001)));
    };

    return (
        <div className="min-h-screen p-4 md:p-6">
            <div className="max-w-6xl mx-auto space-y-6">
                {/* Header */}
                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    className="text-center"
                >
                    <h1 className="text-4xl font-bold mb-2">
                        <span className="neon-text-cyan">👁️ Holo-View</span>{' '}
                        <span className="text-white">AR</span>
                    </h1>
                    <p className="text-slate-400">See molecules floating in your real room</p>
                </motion.div>

                {/* Molecule Selector */}
                <div className="flex flex-wrap justify-center gap-3">
                    {MOLECULES.map((mol) => (
                        <button
                            key={mol.id}
                            onClick={() => { setSelectedMolecule(mol); awardXP(3, `Viewed ${mol.name}`, 'Holo-View'); }}
                            className="scholar-btn px-4 py-2 rounded-xl text-sm cursor-pointer"
                            style={{
                                background: selectedMolecule.id === mol.id ? 'rgba(0, 240, 255, 0.2)' : 'rgba(255,255,255,0.03)',
                                borderColor: selectedMolecule.id === mol.id ? '#00f0ff' : 'rgba(255,255,255,0.1)',
                            }}
                        >
                            🧬 {mol.name.split(' ')[0]}
                        </button>
                    ))}
                </div>

                {/* Viewer */}
                <div className="relative rounded-2xl overflow-hidden" style={{ height: '60vh', border: '2px solid rgba(0, 240, 255, 0.2)' }}>
                    {/* Webcam Background */}
                    <AnimatePresence>
                        {showWebcam && (
                            <motion.div
                                initial={{ opacity: 0 }}
                                animate={{ opacity: 1 }}
                                exit={{ opacity: 0 }}
                                className="absolute inset-0 z-0"
                            >
                                <Webcam
                                    audio={false}
                                    mirrored
                                    className="w-full h-full object-cover"
                                    videoConstraints={{ facingMode: 'user' }}
                                />
                                {/* AR overlay tint */}
                                <div className="absolute inset-0" style={{ background: 'rgba(0, 10, 20, 0.2)' }} />
                            </motion.div>
                        )}
                    </AnimatePresence>

                    {/* No-webcam background */}
                    {!showWebcam && (
                        <div className="absolute inset-0 holo-grid" style={{ background: 'var(--scholar-bg)' }} />
                    )}

                    {/* 3D Canvas */}
                    <div className="absolute inset-0 z-10">
                        <canvas
                            ref={canvasRef}
                            className="w-full h-full cursor-grab active:cursor-grabbing"
                            onMouseDown={handleMouseDown}
                            onMouseMove={handleMouseMove}
                            onMouseUp={() => setIsDragging(false)}
                            onMouseLeave={() => setIsDragging(false)}
                            onWheel={handleWheel}
                        />
                    </div>

                    {/* Info Overlay */}
                    <AnimatePresence>
                        {showInfo && (
                            <motion.div
                                initial={{ opacity: 0, x: -20 }}
                                animate={{ opacity: 1, x: 0 }}
                                exit={{ opacity: 0, x: -20 }}
                                className="absolute bottom-4 left-4 z-20 scholar-card p-4 max-w-sm"
                                style={{ background: 'rgba(10, 14, 26, 0.9)' }}
                            >
                                <h3 className="text-sm font-bold text-white mb-1">{selectedMolecule.name}</h3>
                                <div className="text-xs neon-text-cyan mb-1">{selectedMolecule.formula}</div>
                                <p className="text-xs text-slate-400 mb-2">{selectedMolecule.description}</p>
                                <p className="text-xs text-slate-300">{selectedMolecule.teachingPoint}</p>
                            </motion.div>
                        )}
                    </AnimatePresence>

                    {/* Controls */}
                    <div className="absolute top-4 right-4 z-20 flex flex-col gap-2">
                        <button
                            onClick={() => setShowWebcam(!showWebcam)}
                            className={`scholar-btn px-3 py-2 text-xs rounded-lg cursor-pointer ${showWebcam ? 'scholar-btn-success' : ''}`}
                        >
                            📷 {showWebcam ? 'AR ON' : 'AR OFF'}
                        </button>
                        <button
                            onClick={() => setShowInfo(!showInfo)}
                            className="scholar-btn px-3 py-2 text-xs rounded-lg cursor-pointer"
                        >
                            ℹ️ Info
                        </button>
                        <button
                            onClick={() => setScale(1)}
                            className="scholar-btn px-3 py-2 text-xs rounded-lg cursor-pointer"
                        >
                            🔄 Reset
                        </button>
                    </div>
                </div>

                {/* Controls Legend */}
                <div className="scholar-card p-4 flex flex-wrap gap-6 justify-center text-xs text-slate-500">
                    <span>🖱️ Drag to rotate</span>
                    <span>🔄 Scroll to zoom</span>
                    <span>📷 Toggle AR camera</span>
                    <span>🧬 Click molecules to switch</span>
                </div>
            </div>
        </div>
    );
}
