// Main Application Logic
document.addEventListener('DOMContentLoaded', () => {
    initNavigation();
    initDrugs();
    initReceptors();
    initEffects();
    initSources();
    initModal();
});

// Navigation
function initNavigation() {
    document.querySelectorAll('.nav-links a').forEach(link => {
        link.addEventListener('click', (e) => {
            e.preventDefault();
            const target = link.getAttribute('href').substring(1);
            document.querySelectorAll('.section').forEach(s => s.classList.remove('active'));
            document.querySelectorAll('.nav-links a').forEach(a => a.classList.remove('active'));
            document.getElementById(target).classList.add('active');
            link.classList.add('active');
            
            // Close any open modal when switching tabs
            const modal = document.getElementById('modal');
            if (modal) modal.classList.remove('active');
            if (window.brain3d) window.brain3d.resetHighlights();
        });
    });
}

// Drugs Section
function initDrugs() {
    const grid = document.getElementById('drugs-grid');
    if (!DATA.drugs) return;

    Object.entries(DATA.drugs).forEach(([name, drug]) => {
        const card = createDrugCard(name, drug);
        grid.appendChild(card);
    });

    // Filter buttons
    document.querySelectorAll('.filter-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            const cls = btn.dataset.class;
            document.querySelectorAll('#drugs-grid .card').forEach(card => {
                card.style.display = (cls === 'all' || card.dataset.class === cls) ? 'block' : 'none';
            });
        });
    });
}

function createDrugCard(name, drug) {
    const card = document.createElement('div');
    card.className = 'card';
    card.dataset.class = drug.class || 'other';

    const receptors = Object.keys(drug.receptors || {}).slice(0, 5);
    
    // Find effects for this drug from scraped data
    let drugEffects = [];
    if (typeof EFFECT_DRUG_RELATIONSHIPS !== 'undefined') {
        const drugNameLower = name.toLowerCase();
        const drugNameNormalized = name.replace(/[^a-zA-Z0-9]/g, '').toLowerCase();
        
        Object.keys(EFFECT_DRUG_RELATIONSHIPS).forEach(effectId => {
            const drugs = EFFECT_DRUG_RELATIONSHIPS[effectId];
            const found = drugs.some(d => {
                const dLower = d.toLowerCase();
                const dNormalized = d.replace(/[^a-zA-Z0-9]/g, '').toLowerCase();
                return dLower === drugNameLower || 
                       dNormalized === drugNameNormalized ||
                       dLower.includes(drugNameLower) ||
                       drugNameLower.includes(dLower);
            });
            
            if (found) {
                let effectName = effectId.replace(/_/g, ' ');
                if (typeof EFFECT_METADATA !== 'undefined' && EFFECT_METADATA[effectId]) {
                    effectName = EFFECT_METADATA[effectId].name;
                }
                drugEffects.push(effectName);
            }
        });
    }
    
    const effectsToShow = drugEffects.slice(0, 4);
    const effectsHtml = effectsToShow.length > 0 
        ? `<div class="effect-list" style="margin-top:0.5rem">
             ${effectsToShow.map(e => `<span class="tag effect-tag">${e}</span>`).join('')}
             ${drugEffects.length > 4 ? `<span class="tag more-tag">+${drugEffects.length - 4}</span>` : ''}
           </div>`
        : '';
    
    card.innerHTML = `
        <div class="card-header">
            <span class="card-title">${name}</span>
            <span class="card-badge ${drug.class}">${drug.class}</span>
        </div>
        <div class="card-content">
            <div class="receptor-list">
                ${receptors.map(r => `<span class="tag">${r}</span>`).join('')}
                ${Object.keys(drug.receptors || {}).length > 5 ? '<span class="tag">+more</span>' : ''}
            </div>
            ${effectsHtml}
        </div>
    `;

    card.addEventListener('click', () => showDrugModal(name, drug));
    return card;
}

function showDrugModal(name, drug) {
    const modal = document.getElementById('modal');
    const body = document.getElementById('modal-body');

    const receptorBars = Object.entries(drug.receptors || {}).map(([rec, energy]) => {
        const pct = Math.min(100, Math.abs(energy) / 5);
        return `
            <div class="binding-bar">
                <span class="binding-label">${rec}</span>
                <div class="binding-value">
                    <div class="binding-fill" style="width: ${pct}%"></div>
                </div>
                <span style="margin-left:8px;font-size:0.8rem">${energy} kJ/mol</span>
            </div>
        `;
    }).join('');

    // Find effects for this drug from scraped data
    let drugEffects = [];
    let effectSource = '';
    
    if (typeof EFFECT_DRUG_RELATIONSHIPS !== 'undefined') {
        // Search through all effects to find ones that include this drug
        const drugNameLower = name.toLowerCase();
        const drugNameNormalized = name.replace(/[^a-zA-Z0-9]/g, '').toLowerCase();
        
        Object.keys(EFFECT_DRUG_RELATIONSHIPS).forEach(effectId => {
            const drugs = EFFECT_DRUG_RELATIONSHIPS[effectId];
            const found = drugs.some(d => {
                const dLower = d.toLowerCase();
                const dNormalized = d.replace(/[^a-zA-Z0-9]/g, '').toLowerCase();
                return dLower === drugNameLower || 
                       dNormalized === drugNameNormalized ||
                       dLower.includes(drugNameLower) ||
                       drugNameLower.includes(dLower);
            });
            
            if (found) {
                // Get display name from metadata
                let effectName = effectId.replace(/_/g, ' ');
                if (typeof EFFECT_METADATA !== 'undefined' && EFFECT_METADATA[effectId]) {
                    effectName = EFFECT_METADATA[effectId].name;
                    if (!effectSource && EFFECT_METADATA[effectId].source) {
                        effectSource = EFFECT_METADATA[effectId].source;
                    }
                }
                drugEffects.push({
                    id: effectId,
                    name: effectName
                });
            }
        });
    }

    // Generate effects HTML
    let effectsHtml = '';
    if (drugEffects.length > 0) {
        const effectTags = drugEffects.slice(0, 20).map(e => 
            `<span class="tag effect-tag">${e.name}</span>`
        ).join('');
        const moreCount = drugEffects.length > 20 ? `<span class="tag more-tag">+${drugEffects.length - 20} more</span>` : '';
        
        effectsHtml = `
            <div class="detail-section">
                <h4>Subjective Effects (${drugEffects.length})</h4>
                <p class="source-note">From ${effectSource || 'PsychonautWiki + TripSit'}</p>
                <div class="effect-list">
                    ${effectTags}${moreCount}
                </div>
            </div>
        `;
    }

    // Generate binding animations HTML
    let animationsHtml = '';
    
    if (typeof GIF_INDEX !== 'undefined') {
        // Try to find GIFs - exact match first
        let gifData = GIF_INDEX[name];
        
        // If not found, try variations
        if (!gifData) {
            const variations = [
                name.replace(/ /g, '_'),
                name.replace(/ /g, '-'),
                name.replace(/[^a-zA-Z0-9-]/g, '')
            ];
            for (const v of variations) {
                if (GIF_INDEX[v]) {
                    gifData = GIF_INDEX[v];
                    break;
                }
            }
        }
        
        if (gifData && gifData.length > 0) {
            // Limit to first 15 receptors for UI clarity
            const displayGifs = gifData.slice(0, 15);
            
            const gifItems = displayGifs.map((g, idx) => `
                <div class="gif-item ${idx === 0 ? 'active' : ''}" data-gif="animations/${g.file}" data-receptor="${g.receptor}">
                    <span class="gif-receptor-label">${g.receptor}</span>
                </div>
            `).join('');
            
            animationsHtml = `
                <div class="detail-section">
                    <h4>🎬 Molecular Dynamics (${gifData.length} receptors)</h4>
                    <div class="gif-viewer">
                        <div class="gif-display">
                            <img id="current-gif" src="animations/${displayGifs[0].file}" alt="Binding animation" onerror="this.style.display='none'">
                            <div class="gif-caption">
                                <span id="gif-receptor-name">${displayGifs[0].receptor}</span> receptor binding
                            </div>
                        </div>
                        <div class="gif-selector">
                            ${gifItems}
                            ${gifData.length > 15 ? `<div class="gif-item" style="opacity:0.6">+${gifData.length - 15} more</div>` : ''}
                        </div>
                    </div>
                </div>
            `;
        }
    }

    // Determine source URLs
    const pwikiUrl = `https://psychonautwiki.org/wiki/${encodeURIComponent(name.replace(/ /g, '_'))}`;
    const tripsitUrl = `https://tripsit.me/factsheets/${encodeURIComponent(name.replace(/ /g, '-').toLowerCase())}`;

    body.innerHTML = `
        <h2>${name}</h2>
        <span class="card-badge ${drug.class}" style="display:inline-block;margin-bottom:1rem">${drug.class}</span>
        
        ${animationsHtml}
        
        ${effectsHtml}
        
        <div class="detail-section">
            <h4>Receptor Binding Profile</h4>
            ${receptorBars || '<p>No binding data available</p>'}
        </div>

        <div class="detail-section">
            <h4>Sources</h4>
            <div class="source-citations">
                <div class="source-citation">
                    <span class="source-label">📚</span>
                    <a href="${pwikiUrl}" target="_blank">PsychonautWiki - ${name}</a>
                </div>
                <div class="source-citation">
                    <span class="source-label">📚</span>
                    <a href="${tripsitUrl}" target="_blank">TripSit - ${name}</a>
                </div>
            </div>
        </div>
    `;

    // Add GIF selector event listeners
    document.querySelectorAll('.gif-item').forEach(item => {
        item.addEventListener('click', () => {
            document.querySelectorAll('.gif-item').forEach(i => i.classList.remove('active'));
            item.classList.add('active');
            document.getElementById('current-gif').src = item.dataset.gif;
            document.getElementById('gif-receptor-name').textContent = item.dataset.receptor;
        });
    });

    modal.classList.add('active');
}

// Receptors Section
function initReceptors() {
    const grid = document.getElementById('receptors-grid');
    if (!DATA.receptors) return;

    Object.entries(DATA.receptors).forEach(([name, rec]) => {
        const card = createReceptorCard(name, rec);
        grid.appendChild(card);
    });
}

function createReceptorCard(name, rec) {
    const card = document.createElement('div');
    card.className = 'card';

    const regions = (rec.regions || []).slice(0, 3);
    const effects = (rec.effects || []).slice(0, 4);

    card.innerHTML = `
        <div class="card-header">
            <span class="card-title">${name}</span>
        </div>
        <div class="card-content">
            <p style="margin-bottom:0.5rem">${rec.function || 'Neuroreceptor'}</p>
            <div class="effect-list">
                ${effects.map(e => `<span class="tag">${e.replace(/_/g, ' ')}</span>`).join('')}
            </div>
            <p style="margin-top:0.75rem;font-size:0.8rem;color:var(--text-secondary)">
                Top regions: ${regions.map(r => r.region).join(', ')}
            </p>
        </div>
    `;

    card.addEventListener('click', () => showReceptorModal(name, rec));
    return card;
}

function showReceptorModal(name, rec) {
    const modal = document.getElementById('modal');
    const body = document.getElementById('modal-body');

    // Get source info
    const source = DATA.receptorSources?.[name] || {};
    const pmid = source.primary?.pmid || 'N/A';
    const authors = source.primary?.authors || 'Unknown';
    const year = source.primary?.year || '';
    const journal = source.primary?.journal || '';

    const regionBars = (rec.regions || []).map(r => {
        const pct = Math.min(100, (r.density / 500) * 100);
        return `
            <div class="binding-bar">
                <span class="binding-label" style="width:180px">${r.region.replace(/_/g, ' ')}</span>
                <div class="binding-value">
                    <div class="binding-fill" style="width: ${pct}%"></div>
                </div>
                <span style="margin-left:8px;font-size:0.8rem">${r.density} fmol/mg</span>
            </div>
        `;
    }).join('');

    body.innerHTML = `
        <h2>${name}</h2>
        <p style="color:var(--text-secondary);margin-bottom:1rem">${rec.function || ''}</p>
        
        <div class="detail-section">
            <h4>Brain Distribution</h4>
            ${regionBars || '<p>No distribution data</p>'}
        </div>

        <div class="detail-section">
            <h4>Associated Effects</h4>
            <div class="effect-list">
                ${(rec.effects || []).map(e => `<span class="tag">${e.replace(/_/g, ' ')}</span>`).join('')}
            </div>
        </div>

        <div class="detail-section">
            <h4>Source</h4>
            <div class="source-citation">
                ${authors} (${year}). ${journal}.<br>
                PMID: <a href="https://pubmed.ncbi.nlm.nih.gov/${pmid}" target="_blank">${pmid}</a>
            </div>
        </div>
    `;

    modal.classList.add('active');

    // Highlight in 3D brain
    if (window.brain3d) {
        window.brain3d.highlightByReceptor(rec);
    }
}

// Effects Section
function initEffects() {
    const grid = document.getElementById('effects-grid');
    if (!DATA.effects) return;

    Object.entries(DATA.effects).forEach(([name, eff]) => {
        const card = createEffectCard(name, eff);
        grid.appendChild(card);
    });
}

function createEffectCard(name, eff) {
    const card = document.createElement('div');
    card.className = 'card';

    const desc = DATA.effectDescriptions?.[name] || '';
    const receptors = (eff.receptors || []).slice(0, 4);
    
    // Find drugs associated with this effect
    let effectDrugs = [];
    if (typeof EFFECT_DRUG_RELATIONSHIPS !== 'undefined') {
        const effectDrugsRaw = EFFECT_DRUG_RELATIONSHIPS[name] || [];
        effectDrugs = effectDrugsRaw.slice(0, 4);
    }
    
    const drugsHtml = effectDrugs.length > 0 
        ? `<div class="drug-list" style="margin-top:0.5rem">
             ${effectDrugs.map(d => `<span class="tag drug-tag">${d}</span>`).join('')}
             ${(EFFECT_DRUG_RELATIONSHIPS[name] || []).length > 4 ? `<span class="tag more-tag">+${(EFFECT_DRUG_RELATIONSHIPS[name] || []).length - 4}</span>` : ''}
           </div>`
        : '';

    card.innerHTML = `
        <div class="card-header">
            <span class="card-title">${name.replace(/_/g, ' ')}</span>
        </div>
        <div class="card-content">
            <p style="margin-bottom:0.75rem">${desc}</p>
            <div class="receptor-list">
                ${receptors.map(r => `<span class="tag">${r}</span>`).join('')}
            </div>
            ${drugsHtml}
        </div>
    `;

    card.addEventListener('click', () => showEffectModal(name, eff));
    return card;
}

function showEffectModal(name, eff) {
    const modal = document.getElementById('modal');
    const body = document.getElementById('modal-body');

    const desc = DATA.effectDescriptions?.[name] || '';
    
    // Find drugs associated with this effect
    let effectDrugs = [];
    if (typeof EFFECT_DRUG_RELATIONSHIPS !== 'undefined') {
        effectDrugs = EFFECT_DRUG_RELATIONSHIPS[name] || [];
    }
    
    const drugsHtml = effectDrugs.length > 0
        ? `<div class="detail-section">
            <h4>Associated Drugs (${effectDrugs.length})</h4>
            <div class="effect-list">
                ${effectDrugs.map(d => `<span class="tag drug-tag">${d}</span>`).join('')}
            </div>
           </div>`
        : '';

    body.innerHTML = `
        <h2>${name.replace(/_/g, ' ')}</h2>
        <p style="color:var(--text-secondary);margin-bottom:1rem">${desc}</p>
        
        ${drugsHtml}

        <div class="detail-section">
            <h4>Source</h4>
            <div class="source-citation">
                Effect descriptions: <a href="https://psychonautwiki.org/wiki/Subjective_effect_index" target="_blank">PsychonautWiki Subjective Effect Index</a>
            </div>
        </div>
    `;

    modal.classList.add('active');
}

// Sources Section
function initSources() {
    const drugSources = document.getElementById('drug-sources');
    if (!drugSources) return;

    // Drug effect sources
    drugSources.innerHTML = `
        <div class="source-item">
            <strong>PsychonautWiki</strong>
            <div class="source-meta">
                Community-sourced subjective effect documentation<br>
                <a href="https://psychonautwiki.org" target="_blank">https://psychonautwiki.org</a>
            </div>
        </div>
        <div class="source-item">
            <strong>TripSit</strong>
            <div class="source-meta">
                Drug combination safety data<br>
                <a href="https://tripsit.me" target="_blank">https://tripsit.me</a>
            </div>
        </div>
        <div class="source-item">
            <strong>Erowid Experience Vaults</strong>
            <div class="source-meta">
                First-person experience reports<br>
                <a href="https://erowid.org/experiences/" target="_blank">https://erowid.org</a>
            </div>
        </div>
    `;
}

// Region Select
function initRegionSelect() {
    const select = document.getElementById('region-select');
    
    const regions = [
        'prefrontal_cortex', 'anterior_cingulate', 'posterior_cingulate',
        'claustrum', 'hippocampus', 'amygdala', 'thalamus', 'hypothalamus',
        'nucleus_accumbens', 'striatum', 'visual_cortex', 'dorsal_raphe',
        'locus_coeruleus', 'vta', 'substantia_nigra', 'cerebellum',
        'insula', 'temporal_cortex', 'parietal_cortex'
    ];

    regions.forEach(r => {
        const opt = document.createElement('option');
        opt.value = r;
        opt.textContent = r.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
        select.appendChild(opt);
    });

    select.addEventListener('change', () => {
        if (select.value && window.brain3d) {
            window.brain3d.highlightRegion(select.value);
            showRegionInfo(select.value);
        }
    });
}

window.showRegionInfo = function(region) {
    const info = document.getElementById('region-info');
    
    // Find receptors in this region
    const receptorsHere = [];
    if (DATA.receptors) {
        Object.entries(DATA.receptors).forEach(([name, rec]) => {
            const found = (rec.regions || []).find(r => 
                r.region.toLowerCase().includes(region.toLowerCase()) ||
                region.toLowerCase().includes(r.region.toLowerCase().split('_')[0])
            );
            if (found) {
                receptorsHere.push({ name, density: found.density });
            }
        });
    }

    receptorsHere.sort((a, b) => (b.density || 0) - (a.density || 0));

    info.innerHTML = `
        <h4 style="color:var(--accent-cyan);margin-bottom:0.75rem">
            ${region.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}
        </h4>
        <p style="font-size:0.9rem;margin-bottom:1rem">Receptors present:</p>
        ${receptorsHere.map(r => `
            <div style="display:flex;justify-content:space-between;padding:0.25rem 0;border-bottom:1px solid rgba(255,255,255,0.1)">
                <span>${r.name}</span>
                <span style="color:var(--accent-purple)">${r.density || '?'} fmol/mg</span>
            </div>
        `).join('') || '<p>No data available</p>'}
    `;
};

// Modal
function initModal() {
    const modal = document.getElementById('modal');
    const close = modal.querySelector('.close');

    close.addEventListener('click', () => {
        modal.classList.remove('active');
        if (window.brain3d) window.brain3d.resetHighlights();
    });

    modal.addEventListener('click', (e) => {
        if (e.target === modal) {
            modal.classList.remove('active');
            if (window.brain3d) window.brain3d.resetHighlights();
        }
    });
}
