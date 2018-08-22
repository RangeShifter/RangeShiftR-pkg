object frmGenerateLand: TfrmGenerateLand
  Left = 0
  Top = 0
  HelpContext = 900
  Caption = 'Artificial Landscape'
  ClientHeight = 673
  ClientWidth = 725
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -13
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  Position = poDesigned
  DesignSize = (
    725
    673)
  PixelsPerInch = 96
  TextHeight = 16
  object GridLandLabel: TLabel
    Left = 16
    Top = 302
    Width = 246
    Height = 16
    Caption = 'Landscape parameters for suitable habitat:'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentFont = False
  end
  object Label10: TLabel
    Left = 16
    Top = 119
    Width = 106
    Height = 16
    AutoSize = False
    Caption = 'Landscape type:'
    Color = clWindowText
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentColor = False
    ParentFont = False
    WordWrap = True
  end
  object Label7: TLabel
    Left = 8
    Top = 8
    Width = 391
    Height = 16
    AutoSize = False
    Caption = 'GENERATE ARTIFICIAL LANDSCAPES (1 suitable habitat type)'
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
  object GridMatrixLabel: TLabel
    Left = 375
    Top = 302
    Width = 240
    Height = 16
    Caption = 'Proportions of second matrix habitat type:'
    Visible = False
  end
  object edtXdim: TLabeledEdit
    Left = 108
    Top = 242
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 75
    EditLabel.Height = 16
    EditLabel.Caption = 'X dimension:'
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 6
    Text = '50'
  end
  object edtYdim: TLabeledEdit
    Left = 108
    Top = 272
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 74
    EditLabel.Height = 16
    EditLabel.Caption = 'Y dimension:'
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 7
    Text = '100'
  end
  object edtMaxProp: TLabeledEdit
    Left = 288
    Top = 242
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 88
    EditLabel.Height = 16
    EditLabel.Caption = 'Max. % of hab:'
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 8
    Text = '100.0'
    Visible = False
  end
  object edtMinProp: TLabeledEdit
    Left = 288
    Top = 272
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 85
    EditLabel.Height = 16
    EditLabel.Caption = 'Min. % of hab:'
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    EditLabel.Transparent = True
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 9
    Text = '1.0'
    Visible = False
  end
  object BtnOK: TButton
    Left = 503
    Top = 623
    Width = 89
    Height = 25
    Caption = 'OK'
    Default = True
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ModalResult = 1
    ParentFont = False
    TabOrder = 16
    OnClick = BtnOKClick
  end
  object edtResolution: TLabeledEdit
    Left = 108
    Top = 212
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 89
    EditLabel.Height = 16
    EditLabel.Caption = 'Resolution (m):'
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 5
    Text = '100'
  end
  object edtSeriesNr: TLabeledEdit
    Left = 295
    Top = 146
    Width = 50
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 57
    EditLabel.Height = 16
    EditLabel.Caption = 'Series nr:'
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 3
    Text = '0'
    Visible = False
  end
  object edtNLand: TLabeledEdit
    Left = 295
    Top = 176
    Width = 50
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 66
    EditLabel.Height = 16
    EditLabel.Caption = 'Nr. of reps:'
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 4
    Text = '0'
    Visible = False
  end
  object StringGridLand: TStringGrid
    Left = 16
    Top = 324
    Width = 330
    Height = 105
    ColCount = 2
    RowCount = 2
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
    ParentFont = False
    TabOrder = 10
    ColWidths = (
      64
      64)
    RowHeights = (
      24
      24)
  end
  object StringGridMatrix: TStringGrid
    Left = 375
    Top = 324
    Width = 330
    Height = 105
    ColCount = 2
    RowCount = 2
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
    ParentFont = False
    TabOrder = 11
    Visible = False
    ColWidths = (
      64
      64)
    RowHeights = (
      24
      24)
  end
  object BtnCreateSeries: TButton
    Left = 16
    Top = 489
    Width = 161
    Height = 25
    Caption = 'Create New Series'
    Enabled = False
    TabOrder = 12
    Visible = False
    OnClick = BtnCreateSeriesClick
  end
  object RGLandSeries: TRadioGroup
    Left = 16
    Top = 30
    Width = 329
    Height = 76
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ItemIndex = 0
    Items.Strings = (
      'Generate a landscape at each replicate'
      'Generate && save a series of landscapes')
    ParentFont = False
    TabOrder = 0
    OnClick = RGLandSeriesClick
  end
  object StatusBar1: TStatusBar
    Left = 0
    Top = 654
    Width = 725
    Height = 19
    Panels = <
      item
        Width = 50
      end>
  end
  object RGLandType: TRadioGroup
    Left = 118
    Top = 141
    Width = 99
    Height = 65
    ItemIndex = 0
    Items.Strings = (
      'Discrete'
      'Continuous')
    TabOrder = 2
    OnClick = RGLandTypeClick
  end
  object BtnCancel: TButton
    Left = 616
    Top = 623
    Width = 89
    Height = 25
    Cancel = True
    Caption = 'Cancel'
    ModalResult = 2
    TabOrder = 13
    OnClick = BtnCancelClick
  end
  object RGFract: TRadioGroup
    Left = 16
    Top = 141
    Width = 99
    Height = 65
    ItemIndex = 0
    Items.Strings = (
      'Random'
      'Fractal')
    TabOrder = 1
    OnClick = RGFractClick
  end
  object Memo2: TMemo
    Left = 16
    Top = 520
    Width = 689
    Height = 97
    TabStop = False
    Anchors = [akLeft, akTop, akRight, akBottom]
    BevelKind = bkFlat
    BorderStyle = bsNone
    Color = clInfoBk
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsItalic]
    ParentFont = False
    ReadOnly = True
    TabOrder = 14
    WordWrap = False
  end
  object CBpatch: TCheckBox
    Left = 16
    Top = 435
    Width = 265
    Height = 17
    Caption = 'Generate patch-based landscapes'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 15
    Visible = False
    WordWrap = True
    OnClick = CBpatchClick
  end
  object edtMaxCells: TLabeledEdit
    Left = 143
    Top = 459
    Width = 49
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 123
    EditLabel.Height = 16
    EditLabel.Caption = 'Max. cells per patch: '
    Enabled = False
    LabelPosition = lpLeft
    NumbersOnly = True
    TabOrder = 17
    Text = '100'
    Visible = False
  end
  object CBmatrix: TCheckBox
    Left = 375
    Top = 240
    Width = 288
    Height = 17
    Caption = 'Add second habitat type in the matrix'
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 19
    Visible = False
    OnClick = CBmatrixClick
  end
end
