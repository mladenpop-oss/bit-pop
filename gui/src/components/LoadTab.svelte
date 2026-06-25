<script>
  import { invoke } from '@tauri-apps/api/core';
  import { listen } from '@tauri-apps/api/event';
  import { open } from '@tauri-apps/plugin-dialog';

  let mode = 'map';
  let index1 = '';
  let index2 = '';
  let index3 = '';
  let index4 = '';
  let include2 = false;
  let include3 = false;
  let include4 = false;
  let readsPath = '';
  let outputPath = '';
  let threads = 4;
  let alignMode = 'xor';
  let useTopN = true;
  let topN = 4;
  let consensusTopN = 2;
  let useChunkPct = false;
  let chunkPct = 0.02;
  let status = '';
  let isRunning = false;
  let progress = 0;
  let logLines = [];

  async function selectIndex(slot) {
    const selected = await open({
      title: 'Select Index File',
      filters: [{
        name: 'Bit-Pop Index',
        extensions: ['bitpop']
      }]
    });
    if (selected) {
      if (slot === 1) index1 = selected;
      else if (slot === 2) index2 = selected;
      else if (slot === 3) index3 = selected;
      else if (slot === 4) index4 = selected;
    }
  }

  async function selectReads() {
    const selected = await open({
      title: 'Select Reads File',
      filters: [{
        name: 'FASTQ',
        extensions: ['fastq', 'fq', 'fastq.gz', 'fq.gz']
      }]
    });
    if (selected) {
      readsPath = selected;
      autoGenerateOutput();
    }
  }

  $: if (readsPath) {
    const parts = readsPath.split(/[\\/]/);
    const fileName = parts[parts.length - 1] || '';
    const nameWithoutExt = fileName.replace(/\.(fastq|fq|gz)(\.fastq|\.fq)?$/, '');
    outputPath = readsPath.replace(/[/\\][^/\\]+$/, '') + '\\' + nameWithoutExt + '.sam';
  }

  $: if (mode === 'map') {
    include2 = false;
    include3 = false;
    include4 = false;
  }

  function handleThreadsChange(e) {
    threads = parseInt(e.target.value);
  }

  function handleTopNChange(e) {
    topN = parseInt(e.target.value);
  }

  function handleConsensusTopNChange(e) {
    consensusTopN = parseInt(e.target.value);
  }

  function handleChunkPctChange(e) {
    chunkPct = parseFloat(e.target.value);
  }

  async function runMapping() {
    if (!index1 || !readsPath) {
      status = 'Please select index 1 and reads';
      return;
    }

    isRunning = true;
    status = '';
    progress = 0;
    logLines = [];

    const unlistenStarted = await listen('run-started', (event) => {
      status = event.payload;
    });

    const unlistenProgress = await listen('run-progress', (event) => {
      progress = event.payload;
    });

    const unlistenLog = await listen('run-log', (event) => {
      logLines = [...logLines, event.payload];
    });

    const unlistenFinished = await listen('run-finished', (event) => {
      status = event.payload;
      isRunning = false;
    });

    const unlistenError = await listen('run-error', (event) => {
      status = event.payload;
      isRunning = false;
    });

    try {
      if (mode === 'map') {
        await invoke('run_map', {
          index: index1,
          reads: readsPath,
          output: outputPath,
          alignMode,
          threads,
          useTopN,
          topN,
          useChunkPct,
          chunkPct,
        });
      } else {
        const indexes = [index1];
        if (include2 && index2) indexes.push(index2);
        if (include3 && index3) indexes.push(index3);
        if (include4 && index4) indexes.push(index4);

        await invoke('run_concon', {
          indexes,
          reads: readsPath,
          output: outputPath,
          threads,
          useTopN,
          topN,
          consensusTopN,
          useChunkPct,
          chunkPct,
        });
      }
    } catch (e) {
      status = e;
      isRunning = false;
    }

    unlistenStarted();
    unlistenProgress();
    unlistenLog();
    unlistenFinished();
    unlistenError();
  }

  $: canRun = index1 && readsPath && (!isRunning);
  $: progressPercent = Math.round(progress * 100);
</script>

<div class="loadtab">
  <div class="mode-toggle">
    <button class="mode-btn {mode === 'map' ? 'active' : ''}" on:click={() => mode = 'map'}>Map</button>
    <button class="mode-btn {mode === 'consensus' ? 'active' : ''}" on:click={() => mode = 'consensus'}>Consensus</button>
  </div>

  <div class="field">
    <label>Index 1</label>
    <div class="path-row">
      <span class="path">{index1 || 'No selection'}</span>
      <button on:click={() => selectIndex(1)}>Select</button>
    </div>
  </div>

  {#if mode === 'consensus'}
    <div class="field consensus-index">
      <label>
        <input type="checkbox" bind:checked={include2} />
        Index 2
      </label>
      <div class="path-row">
        <span class="path">{index2 || 'No selection'}</span>
        <button on:click={() => selectIndex(2)} disabled={!include2}>Select</button>
      </div>
    </div>

    <div class="field consensus-index">
      <label>
        <input type="checkbox" bind:checked={include3} />
        Index 3
      </label>
      <div class="path-row">
        <span class="path">{index3 || 'No selection'}</span>
        <button on:click={() => selectIndex(3)} disabled={!include3}>Select</button>
      </div>
    </div>

    <div class="field consensus-index">
      <label>
        <input type="checkbox" bind:checked={include4} />
        Index 4
      </label>
      <div class="path-row">
        <span class="path">{index4 || 'No selection'}</span>
        <button on:click={() => selectIndex(4)} disabled={!include4}>Select</button>
      </div>
    </div>
  {/if}

  <div class="field">
    <label>Reads</label>
    <div class="path-row">
      <span class="path">{readsPath || 'No selection'}</span>
      <button on:click={selectReads}>Select</button>
    </div>
  </div>

  <div class="field">
    <label>Output</label>
    <div class="path-row">
      <input type="text" bind:value={outputPath} class="output-input" placeholder="Auto-generated from reads path" />
    </div>
  </div>

  <div class="field">
    <label>Threads: {threads}</label>
    <div class="slider-row">
      <input type="range" min="1" max="16" value={threads} on:input={handleThreadsChange} />
      <input type="number" min="1" max="16" value={threads} on:input={handleThreadsChange} class="num" />
    </div>
  </div>

  {#if mode === 'map'}
    <div class="field">
      <label>Align Mode</label>
      <select bind:value={alignMode} class="select">
        <option value="xor">xor</option>
        <option value="hybrid">hybrid</option>
        <option value="sw">sw</option>
      </select>
    </div>
  {/if}

  <div class="field">
    <label>
      <input type="checkbox" bind:checked={useTopN} />
      Top-N: {topN}
    </label>
    {#if useTopN}
      <div class="slider-row">
        <input type="range" min="1" max="8" value={topN} on:input={handleTopNChange} />
        <input type="number" min="1" max="8" value={topN} on:input={handleTopNChange} class="num" />
      </div>
    {/if}
  </div>

  {#if mode === 'consensus' && useTopN}
    <div class="field">
      <label>Consensus Top-N: {consensusTopN}</label>
      <div class="slider-row">
        <input type="range" min="2" max="4" value={consensusTopN} on:input={handleConsensusTopNChange} />
        <input type="number" min="2" max="4" value={consensusTopN} on:input={handleConsensusTopNChange} class="num" />
      </div>
    </div>
  {/if}

  <div class="field">
    <label>
      <input type="checkbox" bind:checked={useChunkPct} />
      Chunk Pct: {(chunkPct * 100).toFixed(1)}%
    </label>
    {#if useChunkPct}
      <div class="slider-row">
        <input type="range" min="0.01" max="0.10" step="0.01" value={chunkPct} on:input={handleChunkPctChange} />
        <input type="number" min="0.01" max="0.10" step="0.01" value={chunkPct} on:input={handleChunkPctChange} class="num" />
      </div>
    {/if}
  </div>

  {#if progress > 0}
    <div class="progress-container">
      <div class="progress-bar" style="width: {progressPercent}%"></div>
      <span class="progress-text">{progressPercent}%</span>
    </div>
  {/if}

  <button class="run-btn" on:click={runMapping} disabled={!canRun}>
    {isRunning ? (mode === 'map' ? 'Mapping...' : 'Running...') : (mode === 'map' ? 'Map' : 'Run Consensus')}
  </button>

  {#if status}
    <div class="status">{status}</div>
  {/if}

  {#if logLines.length > 0}
    <div class="log-container">
      {#each logLines as line}
        <div class="log-line">{line}</div>
      {/each}
    </div>
  {/if}
</div>

<style>
  .loadtab {
    display: flex;
    flex-direction: column;
    gap: 15px;
  }
  .mode-toggle {
    display: flex;
    gap: 0;
    margin-bottom: 5px;
    background: #16213e;
    border-radius: 6px;
    padding: 4px;
  }
  .mode-btn {
    flex: 1;
    background: transparent;
    color: #888;
    border: none;
    padding: 8px 16px;
    border-radius: 4px;
    cursor: pointer;
    font-size: 0.9em;
    transition: all 0.2s;
  }
  .mode-btn:hover {
    color: #e0e0e0;
  }
  .mode-btn.active {
    background: #0f3460;
    color: #64ffda;
  }
  .field {
    display: flex;
    flex-direction: column;
    gap: 5px;
  }
  .field label {
    font-size: 0.85em;
    color: #888;
    display: flex;
    align-items: center;
    gap: 8px;
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
  .output-input {
    flex: 1;
    background: #16213e;
    border: 1px solid #333;
    color: #e0e0e0;
    padding: 8px 12px;
    border-radius: 4px;
    font-size: 0.9em;
    font-family: inherit;
  }
  .consensus-index {
    opacity: 0.7;
    padding-left: 20px;
    border-left: 2px solid #333;
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
  .select {
    background: #16213e;
    border: 1px solid #333;
    color: #e0e0e0;
    padding: 8px 12px;
    border-radius: 4px;
    font-size: 0.9em;
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
  .run-btn {
    background: #64ffda;
    color: #1a1a2e;
    font-weight: bold;
    padding: 12px;
    margin-top: 10px;
  }
  .run-btn:hover:not(:disabled) {
    background: #4de8c5;
  }
  .progress-container {
    position: relative;
    height: 24px;
    background: #16213e;
    border-radius: 4px;
    overflow: hidden;
  }
  .progress-bar {
    height: 100%;
    background: #64ffda;
    transition: width 0.3s;
  }
  .progress-text {
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    color: #1a1a2e;
    font-weight: bold;
    font-size: 0.85em;
  }
  .status {
    margin-top: 10px;
    padding: 10px;
    background: #16213e;
    border-radius: 4px;
    font-size: 0.9em;
    text-align: center;
  }
  .log-container {
    margin-top: 10px;
    padding: 10px;
    background: #0a0a1a;
    border-radius: 4px;
    font-size: 0.8em;
    font-family: 'Consolas', 'Courier New', monospace;
    max-height: 200px;
    overflow-y: auto;
  }
  .log-line {
    padding: 2px 0;
    color: #888;
  }
</style>