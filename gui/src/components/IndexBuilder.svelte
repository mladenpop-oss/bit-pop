<script>
   import { invoke } from '@tauri-apps/api/core';
   import { listen } from '@tauri-apps/api/event';
   import { open } from '@tauri-apps/plugin-dialog';

   let genomesPath = '';
   let outputPath = '';
   let k = 10;
   let threads = 4;
   let useCami = false;
   let status = '';
   let isBuilding = false;
   let progress = 0;

   async function selectGenomes() {
     try {
       const selected = await open({
         directory: true,
         multiple: false,
         title: 'Select Genomes Folder'
       });
       if (selected) {
         genomesPath = selected;
       }
     } catch (e) {
       console.error(e);
     }
   }

   $: if (genomesPath) {
     const parts = genomesPath.split(/[\\/]/);
     const folderName = parts[parts.length - 1] || '';
     const shortName = folderName.substring(0, 3);
     const defaultName = `index_${shortName}_k${k}.bitpop`;
     if (!outputPath) {
       outputPath = genomesPath + '\\' + defaultName;
     }
   }
   $: outputPathFixed = outputPath.replace(/k\d+\.bitpop/, `k${k}.bitpop`);
   $: progressPercent = Math.round(progress * 100);

   async function buildIndex() {
     if (!genomesPath || !outputPathFixed) {
       status = 'Please select genomes and output path';
       return;
     }

     isBuilding = true;
     status = '';
     progress = 0;

     listen('run-started', (event) => {
       status = event.payload;
     });

     listen('run-finished', (event) => {
       status = event.payload;
       isBuilding = false;
       progress = 1;
     });

     listen('run-error', (event) => {
       status = event.payload;
       isBuilding = false;
       progress = 0;
     });

     listen('run-progress', (event) => {
       progress = event.payload;
     });

     try {
       await invoke('build_index', {
         genomes: genomesPath,
         output: outputPathFixed,
         k,
         threads,
         cami: useCami,
       });
     } catch (e) {
       status = e;
       isBuilding = false;
       progress = 0;
     }
   }

   function handleKChange(e) {
     k = parseInt(e.target.value);
   }

   function handleThreadsChange(e) {
     threads = parseInt(e.target.value);
   }
</script>

<div class="builder">
  <div class="field">
    <label>Genomes</label>
    <div class="path-row">
      <span class="path">{genomesPath || 'No selection'}</span>
      <button on:click={selectGenomes}>Select</button>
    </div>
  </div>

  <div class="field">
    <label>Output</label>
    <div class="path-row">
      <span class="path">{outputPathFixed || 'Auto-generated'}</span>
    </div>
  </div>

  <div class="field">
    <label>K-mer size: {k}</label>
    <div class="slider-row">
      <input type="range" min="8" max="22" value={k} on:input={handleKChange} />
      <input type="number" min="8" max="22" value={k} on:input={handleKChange} class="num" />
    </div>
  </div>

 <div class="field">
     <label>Threads: {threads}</label>
     <div class="slider-row">
       <input type="range" min="1" max="16" value={threads} on:input={handleThreadsChange} />
       <input type="number" min="1" max="16" value={threads} on:input={handleThreadsChange} class="num" />
     </div>
   </div>

   <div class="field checkbox-field">
     <label class="checkbox-label">
       <input type="checkbox" bind:checked={useCami} />
       <span>--cami</span>
     </label>
   </div>

  <button class="build-btn" on:click={buildIndex} disabled={isBuilding || !genomesPath}>
     {isBuilding ? 'Building...' : 'Build Index'}
   </button>

   {#if isBuilding}
     <div class="progress-container">
       <div class="progress-bar" style="width: {progressPercent}%"></div>
       <span class="progress-text">{progressPercent}%</span>
     </div>
   {/if}

   {#if status}
     <div class="status">{status}</div>
   {/if}
</div>

<style>
  .builder {
    display: flex;
    flex-direction: column;
    gap: 15px;
  }
  .field {
    display: flex;
    flex-direction: column;
    gap: 5px;
  }
  .field label {
    font-size: 0.85em;
    color: #888;
  }
  .path-row {
    display: flex;
    gap: 8px;
    align-items: center;
  }
  .path {
    flex: 1;
    background: #16213e;
    padding: 8px 12px;
    border-radius: 4px;
    font-size: 0.9em;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  }
  .slider-row {
    display: flex;
    gap: 10px;
    align-items: center;
  }
  .slider-row input[type="range"] {
    flex: 1;
  }
  .num {
    width: 50px;
    background: #16213e;
    border: 1px solid #333;
    color: #e0e0e0;
    padding: 4px 8px;
    border-radius: 4px;
    text-align: center;
  }
  button {
    background: #0f3460;
    color: #e0e0e0;
    border: none;
    padding: 8px 16px;
    border-radius: 4px;
    cursor: pointer;
    white-space: nowrap;
  }
  button:disabled {
    opacity: 0.5;
    cursor: not-allowed;
  }
  button:hover:not(:disabled) {
    background: #1a4a7a;
  }
  .build-btn {
    background: #64ffda;
    color: #1a1a2e;
    font-weight: bold;
    padding: 12px;
    margin-top: 10px;
  }
  .build-btn:hover:not(:disabled) {
    background: #4de8c5;
  }
 .status {
     margin-top: 10px;
     padding: 10px;
     background: #16213e;
     border-radius: 4px;
     font-size: 0.9em;
     text-align: center;
   }
   .progress-container {
     position: relative;
     height: 24px;
     background: #16213e;
     border-radius: 4px;
     overflow: hidden;
     margin-top: 10px;
   }
   .progress-bar {
     height: 100%;
     background: #64ffda;
     transition: width 0.3s ease;
   }
   .progress-text {
     position: absolute;
     top: 50%;
     left: 50%;
     transform: translate(-50%, -50%);
     font-size: 0.8em;
     font-weight: bold;
     color: #e0e0e0;
   }
   .checkbox-field {
     display: flex;
     align-items: center;
   }
   .checkbox-label {
     display: flex;
     align-items: center;
     gap: 8px;
     cursor: pointer;
     font-size: 0.9em;
     color: #e0e0e0;
   }
   .checkbox-label input[type="checkbox"] {
     width: 16px;
     height: 16px;
     cursor: pointer;
   }
</style>
