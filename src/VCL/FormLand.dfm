object frmLand: TfrmLand
  Left = 0
  Top = 0
  Caption = 'Landscape'
  ClientHeight = 421
  ClientWidth = 546
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -13
  Font.Name = 'Tahoma'
  Font.Style = [fsBold]
  OldCreateOrder = False
  Position = poDesigned
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 16
  object LabelResolution: TLabel
    Left = 142
    Top = 8
    Width = 89
    Height = 16
    Caption = 'Resolution (m):'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clBlue
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentFont = False
  end
  object LabelNhab: TLabel
    Left = 142
    Top = 60
    Width = 79
    Height = 32
    Caption = 'Nr. of habitat types:'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clBlue
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentFont = False
    Visible = False
    WordWrap = True
  end
  object LabelSpDist: TLabel
    Left = 8
    Top = 300
    Width = 223
    Height = 16
    AutoSize = False
    Caption = 'SPECIES'#39' DISTRIBUTION   '
    Color = clWindowText
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clBlue
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold, fsItalic]
    ParentColor = False
    ParentFont = False
    WordWrap = True
  end
  object Bevel1: TBevel
    Left = -9
    Top = 285
    Width = 570
    Height = 9
    Shape = bsTopLine
  end
  object LabelSpResol: TLabel
    Left = 8
    Top = 328
    Width = 89
    Height = 16
    Caption = 'Resolution (m):'
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clBlue
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentFont = False
  end
  object edtRes: TEdit
    Left = 142
    Top = 30
    Width = 70
    Height = 24
    Alignment = taRightJustify
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentFont = False
    TabOrder = 2
    Text = '0'
  end
  object BtnImport: TButton
    Left = 8
    Top = 254
    Width = 145
    Height = 25
    Caption = 'Import Landscape'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 6
    OnClick = BtnImportClick
  end
  object BtnDynamic: TButton
    Left = 176
    Top = 254
    Width = 169
    Height = 25
    Caption = 'Dynamic Landscape'
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 12
    Visible = False
    OnClick = BtnDynamicClick
  end
  object edtNhabitat: TLabeledEdit
    Left = 142
    Top = 96
    Width = 70
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 4
    EditLabel.Height = 16
    EditLabel.Caption = ' '
    EditLabel.Color = clBlack
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentColor = False
    EditLabel.ParentFont = False
    EditLabel.WordWrap = True
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 3
    Text = '1'
    Visible = False
    OnExit = edtNhabitatExit
  end
  object RGhab: TRadioGroup
    Left = 8
    Top = 8
    Width = 119
    Height = 81
    Caption = 'Raster Type'
    Ctl3D = True
    DoubleBuffered = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ItemIndex = 0
    Items.Strings = (
      'habitat codes'
      '% cover'
      'habitat quality')
    ParentCtl3D = False
    ParentDoubleBuffered = False
    ParentFont = False
    TabOrder = 0
    OnClick = RGhabClick
  end
  object BtnOK: TButton
    Left = 364
    Top = 381
    Width = 75
    Height = 25
    Caption = 'OK'
    TabOrder = 11
    OnClick = BtnOKClick
  end
  object BtnCancel: TButton
    Left = 463
    Top = 381
    Width = 75
    Height = 25
    Cancel = True
    Caption = 'Cancel'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ModalResult = 2
    ParentFont = False
    TabOrder = 10
    OnClick = BtnCancelClick
  end
  object StringGridHab: TStringGrid
    Left = 230
    Top = 30
    Width = 308
    Height = 209
    ColCount = 4
    DefaultColWidth = 60
    DoubleBuffered = False
    Enabled = False
    RowCount = 2
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -11
    Font.Name = 'Tahoma'
    Font.Style = []
    Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
    ParentDoubleBuffered = False
    ParentFont = False
    TabOrder = 4
    ColWidths = (
      60
      60
      60
      60)
    RowHeights = (
      24
      24)
  end
  object RGCellPatch: TRadioGroup
    Left = 8
    Top = 104
    Width = 119
    Height = 72
    Caption = 'Model Type'
    Ctl3D = True
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ItemIndex = 0
    Items.Strings = (
      'Cell-based'
      'Patch-based')
    ParentCtl3D = False
    ParentFont = False
    TabOrder = 1
    WordWrap = True
    OnClick = RGCellPatchClick
  end
  object CBVisualPatch: TCheckBox
    Left = 8
    Top = 182
    Width = 218
    Height = 16
    Caption = 'Visualise patch landscape'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentFont = False
    TabOrder = 5
    Visible = False
    WordWrap = True
  end
  object BtnImportSp: TButton
    Left = 8
    Top = 381
    Width = 216
    Height = 25
    Caption = 'Import Species Distribution '
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 9
    OnClick = BtnImportSpClick
  end
  object edtSpResol: TEdit
    Left = 8
    Top = 351
    Width = 70
    Height = 24
    Alignment = taRightJustify
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentFont = False
    TabOrder = 8
    Text = '0'
  end
  object BtnChangeColours: TButton
    Left = 364
    Top = 254
    Width = 174
    Height = 25
    Caption = 'Change Habitat Colours'
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 7
    Visible = False
    OnClick = BtnChangeColoursClick
  end
  object OpenTextFileDialog1: TOpenTextFileDialog
    Filter = 'Text files(*.txt)|*.txt'
    Left = 272
    Top = 304
  end
end
