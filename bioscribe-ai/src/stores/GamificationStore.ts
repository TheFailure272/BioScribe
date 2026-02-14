'use client';

import { create } from 'zustand';
import { persist } from 'zustand/middleware';

// ─── Badge Definitions ───────────────────────────────────────────────
export interface Badge {
  id: string;
  name: string;
  icon: string;
  description: string;
  unlockedAt?: number; // timestamp
}

export const ALL_BADGES: Badge[] = [
  { id: 'first_drug', name: 'First Synthesis', icon: '🧪', description: 'Design your first drug candidate' },
  { id: 'safe_designer', name: 'Safe Designer', icon: '🛡️', description: 'Pass all Lipinski rules on first try' },
  { id: 'speed_demon', name: 'Speed Demon', icon: '⚡', description: 'Complete Emergency Room under 30 seconds' },
  { id: 'polypharmacy_master', name: 'Polypharmacy Master', icon: '💊', description: 'Identify 5 severe drug interactions' },
  { id: 'organ_protector', name: 'Organ Protector', icon: '🫀', description: 'Achieve zero organ toxicity flags' },
  { id: 'quiz_champion', name: 'Quiz Champion', icon: '🏆', description: 'Score 100% in Discovery Quiz' },
  { id: 'voice_commander', name: 'Voice Commander', icon: '🎙️', description: 'Use 10 Jarvis voice commands' },
  { id: 'molecule_architect', name: 'Molecule Architect', icon: '🧱', description: 'Build 3 stable molecules in Minecraft mode' },
  { id: 'team_player', name: 'Team Player', icon: '🤝', description: 'Collaborate in Metaverse Lab' },
  { id: 'professor_favorite', name: "Professor's Favorite", icon: '🎓', description: 'Read 20 AI Professor explanations' },
  // ─── New Apex Badges ───────────────────────────────────────────────
  { id: 'team_resuscitator', name: 'Team Resuscitator', icon: '🚑', description: 'Complete a Team Code Blue scenario' },
  { id: 'code_blue_leader', name: 'Code Blue Leader', icon: '👨‍⚕️', description: 'Lead as Physician in 3 Code Blue sessions' },
  { id: 'socratic_thinker', name: 'Socratic Thinker', icon: '🏥', description: 'Complete 5 Dr. House cases' },
  { id: 'citizen_scientist', name: 'Citizen Scientist', icon: '🌍', description: 'Design a novel compound submitted for analysis' },
  { id: 'career_ready', name: 'Career Ready', icon: '💼', description: 'Generate your first Career Passport certificate' },
];

// ─── Level Tiers ─────────────────────────────────────────────────────
export interface LevelTier {
  level: number;
  title: string;
  minXP: number;
  icon: string;
  color: string;
}

export const LEVEL_TIERS: LevelTier[] = [
  { level: 0, title: 'Intern', minXP: 0, icon: '🔬', color: '#64748b' },
  { level: 1, title: 'Resident', minXP: 100, icon: '🩺', color: '#3b82f6' },
  { level: 2, title: 'Fellow', minXP: 300, icon: '⚗️', color: '#8b5cf6' },
  { level: 3, title: 'Attending', minXP: 600, icon: '🧬', color: '#f59e0b' },
  { level: 4, title: 'Chief', minXP: 1000, icon: '👨‍⚕️', color: '#ef4444' },
  { level: 5, title: 'Professor', minXP: 2000, icon: '🎓', color: '#ec4899' },
];

// ─── XP Event Log ────────────────────────────────────────────────────
export interface XPEvent {
  amount: number;
  reason: string;
  timestamp: number;
  module: string;
}

// ─── Store Interface ─────────────────────────────────────────────────
interface GamificationState {
  // State
  xp: number;
  level: number;
  badges: Badge[];
  streakDays: number;
  completedChallenges: string[];
  xpHistory: XPEvent[];
  erBestTime: number | null; // seconds
  moleculesBuilt: number;
  voiceCommandsUsed: number;
  explanationsRead: number;
  interactionsFound: number;
  // ─── New Stats ───────────────────────────────────────────────────
  codeBlueCompleted: number;
  houseCasesCompleted: number;
  certificatesGenerated: number;

  // Actions
  awardXP: (amount: number, reason: string, module: string) => void;
  unlockBadge: (badgeId: string) => void;
  addCompletedChallenge: (challengeId: string) => void;
  setERBestTime: (time: number) => void;
  incrementMoleculesBuilt: () => void;
  incrementVoiceCommands: () => void;
  incrementExplanationsRead: () => void;
  incrementInteractionsFound: () => void;
  incrementCodeBlue: () => void;
  incrementHouseCases: () => void;
  incrementCertificates: () => void;
  resetProgress: () => void;
  nukeAll: () => void;
}

// ─── Helper: Calculate Level ─────────────────────────────────────────
function calculateLevel(xp: number): number {
  let level = 0;
  for (const tier of LEVEL_TIERS) {
    if (xp >= tier.minXP) level = tier.level;
  }
  return level;
}

// ─── Zustand Store ───────────────────────────────────────────────────
export const useGamificationStore = create<GamificationState>()(
  persist(
    (set, get) => ({
      // Initial state
      xp: 0,
      level: 0,
      badges: [],
      streakDays: 0,
      completedChallenges: [],
      xpHistory: [],
      erBestTime: null,
      moleculesBuilt: 0,
      voiceCommandsUsed: 0,
      explanationsRead: 0,
      interactionsFound: 0,
      codeBlueCompleted: 0,
      houseCasesCompleted: 0,
      certificatesGenerated: 0,

      // Actions
      awardXP: (amount, reason, module) => {
        set((state) => {
          const newXP = state.xp + amount;
          const newLevel = calculateLevel(newXP);
          return {
            xp: newXP,
            level: newLevel,
            xpHistory: [
              ...state.xpHistory,
              { amount, reason, timestamp: Date.now(), module },
            ],
          };
        });
      },

      unlockBadge: (badgeId) => {
        const state = get();
        if (state.badges.some((b) => b.id === badgeId)) return; // already unlocked
        const badge = ALL_BADGES.find((b) => b.id === badgeId);
        if (!badge) return;
        set({
          badges: [...state.badges, { ...badge, unlockedAt: Date.now() }],
        });
      },

      addCompletedChallenge: (challengeId) => {
        const state = get();
        if (state.completedChallenges.includes(challengeId)) return;
        set({
          completedChallenges: [...state.completedChallenges, challengeId],
        });
      },

      setERBestTime: (time) => {
        const state = get();
        if (state.erBestTime === null || time < state.erBestTime) {
          set({ erBestTime: time });
        }
      },

      incrementMoleculesBuilt: () =>
        set((s) => ({ moleculesBuilt: s.moleculesBuilt + 1 })),

      incrementVoiceCommands: () =>
        set((s) => ({ voiceCommandsUsed: s.voiceCommandsUsed + 1 })),

      incrementExplanationsRead: () =>
        set((s) => ({ explanationsRead: s.explanationsRead + 1 })),

      incrementInteractionsFound: () =>
        set((s) => ({ interactionsFound: s.interactionsFound + 1 })),

      incrementCodeBlue: () =>
        set((s) => ({ codeBlueCompleted: s.codeBlueCompleted + 1 })),

      incrementHouseCases: () =>
        set((s) => ({ houseCasesCompleted: s.houseCasesCompleted + 1 })),

      incrementCertificates: () =>
        set((s) => ({ certificatesGenerated: s.certificatesGenerated + 1 })),

      resetProgress: () =>
        set({
          xp: 0,
          level: 0,
          badges: [],
          streakDays: 0,
          completedChallenges: [],
          xpHistory: [],
          erBestTime: null,
          moleculesBuilt: 0,
          voiceCommandsUsed: 0,
          explanationsRead: 0,
          interactionsFound: 0,
          codeBlueCompleted: 0,
          houseCasesCompleted: 0,
          certificatesGenerated: 0,
        }),

      // ─── Dev Nuke: Clear ALL localStorage ────────────────────────
      nukeAll: () => {
        localStorage.removeItem('bioscribe-scholar-gamification');
        localStorage.removeItem('bioscribe-crowdsourced-cures');
        localStorage.removeItem('bioscribe-metaverse-lab');
        window.location.reload();
      },
    }),
    {
      name: 'bioscribe-scholar-gamification',
    }
  )
);

// ─── Helper Hooks ────────────────────────────────────────────────────
export function getCurrentTier(level: number): LevelTier {
  return LEVEL_TIERS[level] || LEVEL_TIERS[0];
}

export function getNextTier(level: number): LevelTier | null {
  return LEVEL_TIERS[level + 1] || null;
}

export function getXPProgress(xp: number, level: number): number {
  const current = LEVEL_TIERS[level];
  const next = LEVEL_TIERS[level + 1];
  if (!next) return 100; // max level
  return Math.min(100, ((xp - current.minXP) / (next.minXP - current.minXP)) * 100);
}
