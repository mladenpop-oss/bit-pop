<script>
  import { open } from '@tauri-apps/plugin-shell';

  let showStep1 = false;
  let showStep2 = false;
  let showStep3 = false;

  function setStep(step) {
    showStep1 = step === 1;
    showStep2 = step === 2;
    showStep3 = step === 3;
  }

  async function openGitHub() {
    await open('https://github.com/mladenpop-oss/bit-pop');
  }
</script>

<div class="helptab">
  <h2>Bit-Pop DNA Classifier</h2>
  <p class="description">Fast XOR-based genome mapper for DNA read classification.</p>

  <div class="steps">
    <div class="step" class:active={showStep1}>
      <button class="step-header" on:click={() => setStep(showStep1 ? 0 : 1)}>
        <span class="step-num">1</span>
        <span class="step-title">Build Index</span>
        <span class="step-arrow">{showStep1 ? '▾' : '▸'}</span>
      </button>
      {#if showStep1}
        <div class="step-content">
          <ul>
            <li>Select folder that contains reference genomes (FASTA files)</li>
            <li>Pick k-mer size (default 10) - larger k = more precise but slower</li>
            <li>Select number of threads for parallel build</li>
            <li>Use <code>--cami</code> flag for CAMI format genomes</li>
            <li>Click <strong>Build</strong> and wait for the index to be created</li>
          </ul>
        </div>
      {/if}
    </div>

    <div class="step" class:active={showStep2}>
      <button class="step-header" on:click={() => setStep(showStep2 ? 0 : 2)}>
        <span class="step-num">2</span>
        <span class="step-title">Map Reads</span>
        <span class="step-arrow">{showStep2 ? '▾' : '▸'}</span>
      </button>
      {#if showStep2}
        <div class="step-content">
          <p class="mode-info">Two modes:</p>
          <ul>
            <li><strong>Map</strong> - select only one index</li>
            <li><strong>Consensus</strong> - select two or more indexes for consensus classification</li>
          </ul>
          <ul>
            <li>Select reads (FASTQ file)</li>
            <li>Use auto-generated output name or choose your own</li>
            <li>Select number of threads for parallel mapping</li>
            <li>Select align mode: <code>xor</code>, <code>hybrid</code>, or <code>sw</code> (Smith-Waterman)</li>
            <li><strong>Top-N</strong> (default 4) - number of best matches to keep per read</li>
            <li><strong>Consensus Top-N</strong> (default 2) - number of top matches used for calculating consensus score in consensus mode</li>
          </ul>
          <p class="run-info">Click <strong>Map</strong> / <strong>Run Consensus</strong> and watch progress in the log window</p>
        </div>
      {/if}
    </div>

    <div class="step" class:active={showStep3}>
      <button class="step-header" on:click={() => setStep(showStep3 ? 0 : 3)}>
        <span class="step-num">3</span>
        <span class="step-title">Results</span>
        <span class="step-arrow">{showStep3 ? '▾' : '▸'}</span>
      </button>
      {#if showStep3}
        <div class="step-content">
          <ul>
            <li>Load your SAM file and wait for output</li>
            <li>See total reads, mapped/unmapped, and per-genome statistics</li>
            <li>Use filters: search read name, min score, status (mapped/unmapped)</li>
            <li>Click column headers to sort</li>
          </ul>
        </div>
      {/if}
    </div>
  </div>

  <div class="footer">
    <button class="github-btn" on:click={openGitHub}>GitHub Repository</button>
  </div>
</div>

<style>
  .helptab {
    display: flex;
    flex-direction: column;
    gap: 15px;
  }

  .helptab h2 {
    color: #64ffda;
    font-size: 1.1em;
    margin: 0;
  }

  .description {
    color: #aaa;
    margin: 0;
    font-size: 0.9em;
  }

  .steps {
    display: flex;
    flex-direction: column;
    gap: 8px;
  }

  .step {
    border: 1px solid #333;
    border-radius: 6px;
    overflow: hidden;
  }

  .step.active {
    border-color: #64ffda;
  }

  .step-header {
    display: flex;
    align-items: center;
    gap: 10px;
    background: transparent;
    border: none;
    color: #e0e0e0;
    padding: 10px 14px;
    cursor: pointer;
    width: 100%;
    font-size: 0.9em;
  }

  .step-header:hover {
    background: #16213e;
  }

  .step-num {
    background: #0f3460;
    color: #64ffda;
    width: 24px;
    height: 24px;
    border-radius: 50%;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 0.8em;
    font-weight: bold;
  }

  .step.active .step-num {
    background: #64ffda;
    color: #1a1a2e;
  }

  .step-title {
    flex: 1;
  }

  .step-arrow {
    color: #888;
  }

  .step-content {
    padding: 0 14px 14px;
    color: #aaa;
    font-size: 0.85em;
  }

  .step-content ul {
    margin: 8px 0;
    padding-left: 20px;
  }

  .step-content li {
    margin-bottom: 4px;
  }

  .step-content code {
    background: #16213e;
    padding: 2px 6px;
    border-radius: 3px;
    font-size: 0.9em;
  }

  .mode-info, .run-info {
    margin: 8px 0;
    color: #888;
    font-style: italic;
  }

  .footer {
    margin-top: 10px;
    padding-top: 16px;
    border-top: 1px solid #333;
    text-align: center;
  }

  .github-btn {
    background: transparent;
    border: none;
    color: #64ffda;
    text-decoration: none;
    font-size: 0.85em;
    cursor: pointer;
    padding: 4px 8px;
  }

  .github-btn:hover {
    text-decoration: underline;
  }
</style>
