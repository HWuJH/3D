import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# 标题
st.title("SMILES 分子 3D 可视化")

# 输入 SMILES 字符串
smiles = st.text_input("输入 SMILES 字符串", "c1ccccc1")  # 默认苯环

# 选择渲染样式
style = st.selectbox("选择渲染样式", ["stick", "sphere", "stick+sphere"])

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
def render_molecule(mol, style="stick"):
    pdb = Chem.MolToPDBBlock(mol)
    view = py3Dmol.view(width=500, height=400)
    view.addModel(pdb, 'pdb')
    
    # 根据选择的样式渲染
    if style == "stick":
        view.setStyle({'stick': {'colorscheme': 'lightgreyCarbon', 'radius': 0.2}})
    elif style == "sphere":
        view.setStyle({'sphere': {'colorscheme': 'CPK', 'scale': 0.5}})
    elif style == "stick+sphere":
        view.setStyle({'sphere': {'colorscheme': 'CPK', 'scale': 0.3}})
        view.addStyle({'stick': {'colorscheme': 'lightgreyCarbon', 'radius': 0.15}})
    
    view.zoomTo()
    # 将视图转换为 HTML
    html = view._make_html()
    return html

# 显示原子信息
def display_atom_info(mol):
    atom_info = []
    for atom in mol.GetAtoms():
        atom_info.append(f"原子 {atom.GetSymbol()} (ID: {atom.GetIdx()})")
    return atom_info

# 主逻辑
if smiles:
    mol = generate_3d_structure(smiles)
    if mol:
        st.write("### 分子 3D 结构")
        html = render_molecule(mol, style)
        # 在 Streamlit 中显示 HTML
        st.components.v1.html(html, width=500, height=400, scrolling=True)
        
        # 显示分子式
        st.write(f"**分子式**: {Chem.MolToSmiles(mol)}")
        
        # 显示原子信息
        st.write("### 原子信息")
        atom_info = display_atom_info(mol)
        for info in atom_info:
            st.write(info)
