<script>
  import IndexBuilder from './components/IndexBuilder.svelte';
  import LoadTab from './components/LoadTab.svelte';
  import ResultsTab from './components/ResultsTab.svelte';
  import HelpTab from './components/HelpTab.svelte';

  let activeTab = 'build';
</script>

<div class="app">
  <h1 class="title">Bit-Pop</h1>

  <div class="tabs">
    <button class="tab {activeTab === 'build' ? 'active' : ''}" on:click={() => activeTab = 'build'}>Build</button>
    <button class="tab {activeTab === 'load' ? 'active' : ''}" on:click={() => activeTab = 'load'}>Load/Run</button>
    <button class="tab {activeTab === 'results' ? 'active' : ''}" on:click={() => activeTab = 'results'}>Results</button>
    <button class="tab {activeTab === 'help' ? 'active' : ''}" on:click={() => activeTab = 'help'}>Help</button>
  </div>

  {#if activeTab === 'build'}
    <IndexBuilder />
  {:else if activeTab === 'load'}
    <LoadTab />
  {:else if activeTab === 'results'}
    <ResultsTab />
  {:else if activeTab === 'help'}
    <HelpTab />
  {/if}
</div>

<style>
  .app {
    padding: 20px;
    max-width: 700px;
    margin: 0 auto;
  }
  .title {
    font-size: 1.2em;
    margin-bottom: 15px;
    color: #64ffda;
    text-align: center;
  }
  .tabs {
    display: flex;
    gap: 0;
    margin-bottom: 20px;
    border-bottom: 2px solid #333;
  }
  .tab {
    background: transparent;
    color: #888;
    border: none;
    padding: 10px 20px;
    cursor: pointer;
    font-size: 0.95em;
    border-bottom: 2px solid transparent;
    margin-bottom: -2px;
    transition: all 0.2s;
  }
  .tab:hover {
    color: #e0e0e0;
  }
  .tab.active {
    color: #64ffda;
    border-bottom: 2px solid #64ffda;
  }
</style>