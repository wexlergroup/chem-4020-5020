import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { viteSingleFile } from 'vite-plugin-singlefile'

export default defineConfig({
  plugins: [react(), viteSingleFile()],
  // The shared studio-kit lies outside this folder; its imports resolve to the packages installed here.
  resolve: { dedupe: ['react', 'react-dom'] },
  server: { fs: { allow: ['..'] } },
  build: {
    outDir: 'dist',
    emptyOutDir: true,
  }
})
