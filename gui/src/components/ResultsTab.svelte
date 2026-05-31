<script>
  import { invoke } from '@tauri-apps/api/core';
  import { open } from '@tauri-apps/plugin-dialog';

  let samPath = '';
  let isLoading = false;
  let status = '';
  let totalReads = 0;
  let mappedReads = 0;
  let unmappedReads = 0;
  let mappedPct = 0;
  let genomeStats = [];
  let detailRows = [];
  let filterText = '';
  let filterMinScore = 0;
  let filterMapped = 'all';
  let sortCol = 'genome_name';
  let sortAsc = true;
  let page = 0;
  const pageSize = 100;

  async function selectSam() {
    const selected = await open({
      title: 'Select SAM File',
      filters: [{
        name: 'SAM',
        extensions: ['sam']
      }]
    });
    if (selected) {
      samPath = selected;
      loadSam();
    }
  }

  async function loadSam() {
    if (!samPath) return;
    isLoading = true;
    status = 'Loading SAM file...';
    genomeStats = [];
    detailRows = [];
    page = 0;

    try {
      const stats = await invoke('parse_sam_stats', { path: samPath });
      totalReads = stats.total;
      mappedReads = stats.mapped;
      unmappedReads = stats.unmapped;
      mappedPct = totalReads > 0 ? Math.round((mappedReads / totalReads) * 100) : 0;
      genomeStats = stats.genomes || [];
      status = 'Loaded!';
      loadPage();
    } catch (e) {
      status = `Error: ${e}`;
    }

    isLoading = false;
  }

  async function loadPage() {
    if (!samPath) return;
    try {
      const result = await invoke('parse_sam_rows', {
        path: samPath,
        page,
        pageSize,
        filterText,
        filterMinScore,
        filterMapped,
        sortCol,
        sortAsc,
      });
      detailRows = result.rows;
    } catch (e) {
      status = `Error loading rows: ${e}`;
    }
  }

  $: filteredGenomes = genomeStats;

  function handleSort(col) {
    if (sortCol === col) {
      sortAsc = !sortAsc;
    } else {
      sortCol = col;
      sortAsc = true;
    }
    loadPage();
  }

  $: canFilter = filterText || filterMinScore > 0 || filterMapped !== 'all';
</script>

<div class="resultstab">
  <div class="field">
    <label>SAM File</label>
    <div class="path-row">
      <span class="path">{samPath || 'No selection'}</span>
      <button on:click={selectSam}>Select</button>
    </div>
  </div>

  {#if totalReads > 0}
    <div class="summary">
      <div class="summary-card">
        <div class="summary-value">{totalReads.toLocaleString()}</div>
        <div class="summary-label">Total Reads</div>
      </div>
      <div class="summary-card mapped">
        <div class="summary-value">{mappedPct}%</div>
        <div class="summary-label">Mapped ({mappedReads.toLocaleString()})</div>
      </div>
      <div class="summary-card unmapped">
        <div class="summary-value">{100 - mappedPct}%</div>
        <div class="summary-label">Unmapped ({unmappedReads.toLocaleString()})</div>
      </div>
    </div>

    {#if genomeStats.length > 0}
      <div class="section">
        <h3>Per-Genome Summary</h3>
        <table class="genome-table">
          <thead>
            <tr>
              <th>Genome</th>
              <th>Mapped</th>
              <th>%</th>
              <th>Avg Score</th>
            </tr>
          </thead>
          <tbody>
            {#each genomeStats as g}
              <tr>
                <td>{g.genome_name}</td>
                <td>{g.mapped.toLocaleString()}</td>
                <td>{g.pct.toFixed(1)}%</td>
                <td>{g.avg_score.toFixed(4)}</td>
              </tr>
            {/each}
          </tbody>
        </table>
      </div>
    {/if}

    <div class="section">
      <h3>Details</h3>
      <div class="filters">
        <input type="text" placeholder="Search read name..." bind:value={filterText} on:input={loadPage} class="filter-input" />
        <div class="filter-group">
          <label>Min Score:</label>
          <input type="number" min="0" max="100" value={filterMinScore} on:input={(e) => filterMinScore = parseFloat(e.target.value) || 0} class="num" />
        </div>
        <div class="filter-group">
          <label>Status:</label>
          <select bind:value={filterMapped} on:change={loadPage} class="select">
            <option value="all">All</option>
            <option value="mapped">Mapped</option>
            <option value="unmapped">Unmapped</option>
          </select>
        </div>
      </div>

      <table class="detail-table">
        <thead>
          <tr>
            <th on:click={() => handleSort('read_name')}>Read Name {sortCol === 'read_name' ? (sortAsc ? '↑' : '↓') : ''}</th>
            <th on:click={() => handleSort('genome_name')}>Genome {sortCol === 'genome_name' ? (sortAsc ? '↑' : '↓') : ''}</th>
            <th on:click={() => handleSort('score')}>Score {sortCol === 'score' ? (sortAsc ? '↑' : '↓') : ''}</th>
            <th>Position</th>
            <th>CIGAR</th>
            <th on:click={() => handleSort('status')}>Status {sortCol === 'status' ? (sortAsc ? '↑' : '↓') : ''}</th>
          </tr>
        </thead>
        <tbody>
          {#each detailRows as r}
            <tr class={r.status === 'unmapped' ? 'unmapped' : ''}>
              <td class="read-name">{r.read_name}</td>
              <td>{r.genome_name}</td>
              <td>{r.score.toFixed(4)}</td>
              <td>{r.position}</td>
              <td class="cigar">{r.cigar}</td>
              <td>{r.status}</td>
            </tr>
          {/each}
        </tbody>
      </table>

      <div class="pagination">
        <button on:click={() => { page--; loadPage(); }} disabled={page === 0}>Previous</button>
        <span>Page {page + 1}</span>
        <button on:click={() => { page++; loadPage(); }} disabled={detailRows.length < pageSize}>Next</button>
      </div>
    </div>
  {/if}

  {#if status}
    <div class="status">{status}</div>
  {/if}
</div>

<style>
  .resultstab {
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
  .summary {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 10px;
  }
  .summary-card {
    background: #16213e;
    padding: 15px;
    border-radius: 6px;
    text-align: center;
  }
  .summary-card.mapped {
    border-top: 3px solid #64ffda;
  }
  .summary-card.unmapped {
    border-top: 3px solid #ff6b6b;
  }
  .summary-value {
    font-size: 1.5em;
    font-weight: bold;
    color: #e0e0e0;
  }
  .summary-label {
    font-size: 0.8em;
    color: #888;
    margin-top: 4px;
  }
  .section {
    margin-top: 10px;
  }
  .section h3 {
    color: #64ffda;
    font-size: 1em;
    margin-bottom: 10px;
  }
  .genome-table, .detail-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.85em;
  }
  .genome-table th, .genome-table td,
  .detail-table th, .detail-table td {
    padding: 6px 8px;
    text-align: left;
    border-bottom: 1px solid #222;
  }
  .genome-table th, .detail-table th {
    color: #888;
    font-weight: normal;
    cursor: pointer;
    user-select: none;
  }
  .genome-table th:hover, .detail-table th:hover {
    color: #e0e0e0;
  }
  .genome-table td, .detail-table td {
    color: #ccc;
  }
  .detail-table tr.unmapped td {
    color: #666;
  }
  .read-name {
    font-family: 'Consolas', 'Courier New', monospace;
    font-size: 0.9em;
  }
  .cigar {
    font-family: 'Consolas', 'Courier New', monospace;
    font-size: 0.85em;
    color: #888;
  }
  .filters {
    display: flex;
    gap: 15px;
    align-items: center;
    margin-bottom: 10px;
  }
  .filter-input {
    flex: 1;
    background: #16213e;
    border: 1px solid #333;
    color: #e0e0e0;
    padding: 6px 10px;
    border-radius: 4px;
    font-size: 0.85em;
  }
  .filter-group {
    display: flex;
    align-items: center;
    gap: 6px;
  }
  .filter-group label {
    font-size: 0.8em;
    color: #888;
  }
  .num {
    width: 60px;
    background: #16213e;
    border: 1px solid #333;
    color: #e0e0e0;
    padding: 4px 8px;
    border-radius: 4px;
    text-align: center;
    font-size: 0.85em;
  }
  .select {
    background: #16213e;
    border: 1px solid #333;
    color: #e0e0e0;
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 0.85em;
  }
  .pagination {
    display: flex;
    justify-content: center;
    align-items: center;
    gap: 15px;
    margin-top: 10px;
  }
  .pagination span {
    color: #888;
    font-size: 0.85em;
  }
  button {
    background: #0f3460;
    color: #e0e0e0;
    border: none;
    padding: 8px 16px;
    border-radius: 4px;
    cursor: pointer;
    white-space: nowrap;
    font-size: 0.85em;
  }
  button:disabled {
    opacity: 0.5;
    cursor: not-allowed;
  }
  button:hover:not(:disabled) {
    background: #1a4a7a;
  }
  .status {
    margin-top: 10px;
    padding: 10px;
    background: #16213e;
    border-radius: 4px;
    font-size: 0.9em;
    text-align: center;
  }
</style>
