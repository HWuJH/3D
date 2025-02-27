import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# 标题
st.title("SMILES 分子 3D 可视化")

# 输入 SMILES 字符串
smiles = st.text_input("输入 SMILES 字符串", "CCO")

# 生成分子 3D 结构
def generate_3d_structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        mol = Chem.AddHs(mol)  # 添加氢原子
        AllChem.EmbedMolecule(mol)  # 生成 3D 坐标
        AllChem.MMFFOptimizeMolecule(mol)  # 优化结构
        return mol
    else:
        st.error("无效的 SMILES 字符串")
        return None

# 渲染分子
def render_molecule(mol):
    pdb = Chem.MolToPDBBlock(mol)
    view = py3Dmol.view(width=400, height=300)
    view.addModel(pdb, 'pdb')
    view.setStyle({'stick': {}})
    view.zoomTo()
    # 将视图转换为 HTML
    html = view._make_html()
    return html

# 主逻辑
if smiles:
    mol = generate_3d_structure(smiles)
    if mol:
        st.write("### 分子 3D 结构")
        html = render_molecule(mol)
        # 在 Streamlit 中显示 HTML
        st.components.v1.html(html, width=400, height=300, scrolling=True)
