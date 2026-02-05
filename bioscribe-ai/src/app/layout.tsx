import type { Metadata } from "next";
import { Geist, Geist_Mono } from "next/font/google";
import Script from "next/script";
import "./globals.css";
import { MainNavigation } from "@/components/MainNavigation";
import { CollaborationProvider } from "@/contexts/CollaborationContext";

const geistSans = Geist({
  variable: "--font-geist-sans",
  subsets: ["latin"],
});

const geistMono = Geist_Mono({
  variable: "--font-geist-mono",
  subsets: ["latin"],
});

export const metadata: Metadata = {
  title: "BioScribe AI - Drug Discovery Platform",
  description: "AI-powered drug discovery and molecular design platform",
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html lang="en">
      <body
        className={`${geistSans.variable} ${geistMono.variable} antialiased`}
      >
        {/* Load 3DMol.js before interactive components */}
        <Script
          src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"
          strategy="beforeInteractive"
        />
        <CollaborationProvider>
          <MainNavigation />
          {children}
        </CollaborationProvider>
      </body>
    </html>
  );
}

