object frmSeeding: TfrmSeeding
  Left = 571
  Top = 0
  Width = 858
  Height = 600
  HorzScrollBar.Range = 798
  VertScrollBar.Range = 500
  Caption = 'Initialisation Rules'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  Menu = MainMenu_Seeding
  OldCreateOrder = False
  Position = poDesigned
  PixelsPerInch = 96
  TextHeight = 13
  object PanelSpDist: TPanel
    Left = 271
    Top = 160
    Width = 537
    Height = 177
    Alignment = taLeftJustify
    BevelInner = bvLowered
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 4
    VerticalAlignment = taAlignTop
    Visible = False
    object LabelListBox: TLabel
      Left = 12
      Top = 71
      Width = 102
      Height = 80
      Caption = 
        'Species presence cells (x && y, at the resolution of the species' +
        ' distribution data):'
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      Visible = False
      WordWrap = True
    end
    object LabelSpDistn: TLabel
      Left = 5
      Top = 2
      Width = 198
      Height = 16
      Caption = 'FROM SPECIES'#39' DISTRIBUTION   '
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clBlue
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold, fsItalic]
      ParentFont = False
    end
    object LBSpPres: TListBox
      Left = 120
      Top = 74
      Width = 100
      Height = 84
      Enabled = False
      ExtendedSelect = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      MultiSelect = True
      ParentFont = False
      TabOrder = 1
      Visible = False
    end
    object edtNSpCells: TLabeledEdit
      Left = 107
      Top = 32
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 92
      EditLabel.Height = 32
      EditLabel.Caption = 'No. of sp. distibution cells:'
      EditLabel.Font.Charset = DEFAULT_CHARSET
      EditLabel.Font.Color = clWindowText
      EditLabel.Font.Height = -13
      EditLabel.Font.Name = 'Tahoma'
      EditLabel.Font.Style = []
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
      TabOrder = 0
      Text = '0'
    end
  end
  object RG1: TRadioGroup
    Left = 8
    Top = 8
    Width = 245
    Height = 129
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clBlack
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    Items.Strings = (
      'Free initialisation'
      'From species'#39' distribution'
      'From initial individuals file'
      'From initialisation file')
    ParentFont = False
    TabOrder = 0
    OnClick = RG1Click
  end
  object RGRules: TRadioGroup
    Left = 8
    Top = 154
    Width = 245
    Height = 112
    Caption = 'Initialise:'
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clBlack
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ItemIndex = 1
    Items.Strings = (
      'SET IN CODE'
      'SET IN CODE'
      'SET IN CODE')
    ParentFont = False
    TabOrder = 2
    OnClick = RGRulesClick
  end
  object RGCells: TRadioGroup
    Left = 8
    Top = 152
    Width = 245
    Height = 177
    Caption = 'Initialise:'
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ItemIndex = 0
    Items.Strings = (
      'SET IN CODE'
      'SET IN CODE'
      'SET IN CODE')
    ParentFont = False
    TabOrder = 1
    Visible = False
    WordWrap = True
    OnClick = RGCellsClick
  end
  object BtnOK: TBitBtn
    Left = 627
    Top = 495
    Width = 88
    Height = 25
    Caption = 'OK'
    Default = True
    ModalResult = 1
    TabOrder = 7
    OnClick = BtnOKClick
  end
  object BtnCancel: TBitBtn
    Left = 721
    Top = 495
    Width = 88
    Height = 25
    Cancel = True
    Caption = 'Cancel'
    ModalResult = 2
    TabOrder = 8
    OnClick = BtnCancelClick
  end
  object BtnSaveInitial: TBitBtn
    Left = 8
    Top = 335
    Width = 197
    Height = 31
    Caption = 'Save Initialisation File '
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    Glyph.Data = {
      36030000424D3603000000000000360000002800000010000000100000000100
      18000000000000030000120B0000120B00000000000000000000FF00FFFF00FF
      FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00
      FFFF00FFFF00FFFF00FFFF00FFFF00FF97433F97433FB59A9BB59A9BB59A9BB5
      9A9BB59A9BB59A9BB59A9B93303097433FFF00FFFF00FFFF00FFFF00FF97433F
      D66868C66060E5DEDF92292A92292AE4E7E7E0E3E6D9DFE0CCC9CC8F201FAF46
      4697433FFF00FFFF00FFFF00FF97433FD06566C25F5FE9E2E292292A92292AE2
      E1E3E2E6E8DDE2E4CFCCCF8F2222AD464697433FFF00FFFF00FFFF00FF97433F
      D06565C15D5DECE4E492292A92292ADFDDDFE1E6E8E0E5E7D3D0D28A1E1EAB44
      4497433FFF00FFFF00FFFF00FF97433FD06565C15B5CEFE6E6EDE5E5E5DEDFE0
      DDDFDFE0E2E0E1E3D6D0D2962A2AB24A4A97433FFF00FFFF00FFFF00FF97433F
      CD6263C86060C96767CC7272CA7271C66969C46464CC6D6CCA6667C55D5DCD65
      6597433FFF00FFFF00FFFF00FF97433FB65553C27B78D39D9CD7A7A5D8A7A6D8
      A6A5D7A09FD5A09FD7A9A7D8ABABCC666797433FFF00FFFF00FFFF00FF97433F
      CC6667F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9CC66
      6797433FFF00FFFF00FFFF00FF97433FCC6667F9F9F9F9F9F9F9F9F9F9F9F9F9
      F9F9F9F9F9F9F9F9F9F9F9F9F9F9CC666797433FFF00FFFF00FFFF00FF97433F
      CC6667F9F9F9CDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDF9F9F9CC66
      6797433FFF00FFFF00FFFF00FF97433FCC6667F9F9F9F9F9F9F9F9F9F9F9F9F9
      F9F9F9F9F9F9F9F9F9F9F9F9F9F9CC666797433FFF00FFFF00FFFF00FF97433F
      CC6667F9F9F9CDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDCDF9F9F9CC66
      6797433FFF00FFFF00FFFF00FF97433FCC6667F9F9F9F9F9F9F9F9F9F9F9F9F9
      F9F9F9F9F9F9F9F9F9F9F9F9F9F9CC666797433FFF00FFFF00FFFF00FFFF00FF
      97433FF9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F9F99743
      3FFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF
      00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF}
    ParentFont = False
    TabOrder = 6
    OnClick = BtnSaveInitialClick
  end
  object StatusBar1: TStatusBar
    Left = 0
    Top = 543
    Width = 842
    Height = 19
    Panels = <
      item
        Width = 50
      end>
    ExplicitTop = 523
  end
  object PanelNInds: TPanel
    Left = 271
    Top = 8
    Width = 538
    Height = 147
    Alignment = taLeftJustify
    BevelInner = bvLowered
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 3
    VerticalAlignment = taAlignTop
    object LabelNinds: TLabel
      Left = 12
      Top = 7
      Width = 134
      Height = 16
      Caption = 'No. of individuals per'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold]
      ParentFont = False
    end
    object LabelPropn: TLabel
      Left = 229
      Top = 7
      Width = 274
      Height = 16
      Caption = 'Proportion of individuals per stage-class   '
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold]
      ParentFont = False
      Visible = False
    end
    object Bevel6: TBevel
      Left = 210
      Top = 7
      Width = 13
      Height = 132
      Shape = bsLeftLine
    end
    object init_unit: TLabel
      Left = 150
      Top = 7
      Width = 21
      Height = 16
      Caption = 'cell'
    end
    object cboNinds: TComboBox
      Left = 12
      Top = 29
      Width = 84
      Height = 24
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ItemIndex = 1
      ParentFont = False
      TabOrder = 0
      Text = 'at half K'
      OnChange = NindsChange
      Items.Strings = (
        'at K'
        'at half K'
        'set value')
    end
    object edtNinds: TEdit
      Left = 102
      Top = 29
      Width = 55
      Height = 24
      Alignment = taRightJustify
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 1
      Text = '0'
      Visible = False
    end
    object RGinitAges: TRadioGroup
      Left = 13
      Top = 56
      Width = 185
      Height = 85
      Caption = 'Initial ages:'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ItemIndex = 2
      Items.Strings = (
        'Lowest age'
        'Randomise'
        'Quasi-equilibrium')
      ParentFont = False
      TabOrder = 3
    end
    object SGinitialStage: TStringGrid
      Left = 229
      Top = 29
      Width = 227
      Height = 72
      BiDiMode = bdLeftToRight
      ColCount = 3
      DefaultColWidth = 50
      DoubleBuffered = False
      Enabled = False
      RowCount = 2
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
      ParentBiDiMode = False
      ParentDoubleBuffered = False
      ParentFont = False
      TabOrder = 2
      Visible = False
      ColWidths = (
        50
        50
        50)
      RowHeights = (
        24
        24)
    end
  end
  object PanelPlus: TPanel
    Left = 271
    Top = 341
    Width = 538
    Height = 148
    BevelInner = bvLowered
    Enabled = False
    TabOrder = 5
    Visible = False
    DesignSize = (
      538
      148)
    object LabelLandRes: TLabel
      Left = 12
      Top = 69
      Width = 106
      Height = 39
      Anchors = []
      AutoSize = False
      Caption = 'at the landscape resolution'
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsItalic]
      ParentFont = False
      WordWrap = True
      ExplicitTop = 68
    end
    object Bevel3: TBevel
      Left = 265
      Top = 11
      Width = 18
      Height = 131
      Anchors = []
      Shape = bsLeftLine
      ExplicitTop = 10
    end
    object CBAdditional: TCheckBox
      Left = 7
      Top = 0
      Width = 241
      Height = 32
      Anchors = []
      Caption = 'SET IN CODE'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold]
      ParentFont = False
      TabOrder = 0
      WordWrap = True
      OnClick = CBAdditionalClick
    end
    object edtAddX: TLabeledEdit
      Left = 38
      Top = 39
      Width = 50
      Height = 24
      Alignment = taRightJustify
      Anchors = []
      AutoSize = False
      EditLabel.Width = 17
      EditLabel.Height = 16
      EditLabel.Caption = 'X: '
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
      TabOrder = 1
      Text = '0'
    end
    object edtAddY: TLabeledEdit
      Left = 120
      Top = 39
      Width = 50
      Height = 24
      Alignment = taRightJustify
      Anchors = []
      EditLabel.Width = 16
      EditLabel.Height = 16
      EditLabel.Caption = 'Y: '
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
      TabOrder = 2
      Text = '0'
    end
    object BtnAdd: TButton
      Left = 180
      Top = 39
      Width = 43
      Height = 25
      Anchors = []
      Caption = 'Add'
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 4
      OnClick = BtnAddClick
    end
    object InitialList: TMemo
      Left = 289
      Top = 13
      Width = 96
      Height = 121
      Anchors = []
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      ReadOnly = True
      ScrollBars = ssVertical
      TabOrder = 5
    end
    object edtPatchID: TLabeledEdit
      Left = 120
      Top = 39
      Width = 50
      Height = 24
      Alignment = taRightJustify
      Anchors = []
      AutoSize = False
      EditLabel.Width = 56
      EditLabel.Height = 16
      EditLabel.Caption = 'Patch ID: '
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
  end
  object PanelFree: TPanel
    Left = 271
    Top = 161
    Width = 538
    Height = 174
    Alignment = taLeftJustify
    BevelInner = bvLowered
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 9
    VerticalAlignment = taAlignTop
    Visible = False
    object LabInitPos: TLabel
      Left = 38
      Top = 24
      Width = 91
      Height = 16
      Caption = 'Initial Position'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold]
      ParentFont = False
    end
    object LabelFreeInit: TLabel
      Left = 7
      Top = 2
      Width = 144
      Height = 16
      Caption = 'FREE INITIALISATION   '
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clBlue
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold, fsItalic]
      ParentFont = False
    end
    object Bevel5: TBevel
      Left = 275
      Top = 10
      Width = 5
      Height = 155
      Shape = bsLeftLine
    end
    object edtMaxY: TLabeledEdit
      Left = 210
      Top = 76
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 43
      EditLabel.Height = 16
      EditLabel.Caption = 'Max. Y:'
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
      TabOrder = 3
      Text = '0'
    end
    object edtMaxX: TLabeledEdit
      Left = 210
      Top = 46
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 44
      EditLabel.Height = 16
      EditLabel.Caption = 'Max. X:'
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
      TabOrder = 1
      Text = '0'
    end
    object edtMinY: TLabeledEdit
      Left = 94
      Top = 76
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 40
      EditLabel.Height = 16
      EditLabel.Caption = 'Min. Y:'
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
      TabOrder = 2
      Text = '0'
    end
    object edtMinX: TLabeledEdit
      Left = 94
      Top = 46
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 41
      EditLabel.Height = 16
      EditLabel.Caption = 'Min. X:'
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
      TabOrder = 0
      Text = '0'
    end
    object edtTotRandomCells: TLabeledEdit
      Left = 210
      Top = 114
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 75
      EditLabel.Height = 16
      EditLabel.Caption = 'SET IN CODE'
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
    object edtInitFreezeYear: TLabeledEdit
      Left = 470
      Top = 13
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 176
      EditLabel.Height = 16
      EditLabel.Caption = 'Freeze initial range until year: '
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
      Text = '0'
    end
    object CBrestrictRange: TCheckBox
      Left = 312
      Top = 48
      Width = 205
      Height = 17
      Alignment = taLeftJustify
      Caption = 'Restrict range to northern front: '
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 8
      OnClick = CBrestrictRangeClick
    end
    object edtRestrictRows: TLabeledEdit
      Left = 470
      Top = 75
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 75
      EditLabel.Height = 16
      EditLabel.Caption = 'No. of rows: '
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
      TabOrder = 6
      Text = '100'
      Visible = False
    end
    object edtRestrictFreq: TLabeledEdit
      Left = 470
      Top = 105
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 113
      EditLabel.Height = 16
      EditLabel.Caption = 'Frequency (years): '
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
      TabOrder = 7
      Text = '10'
      Visible = False
    end
    object edtFinalFreezeYear: TLabeledEdit
      Left = 470
      Top = 135
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 145
      EditLabel.Height = 16
      EditLabel.Caption = 'Freeze range after year: '
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
      TabOrder = 9
      Text = '1000'
      Visible = False
    end
  end
  object OpenInitialisationFile: TOpenTextFileDialog
    Left = 208
    Top = 373
  end
  object MainMenu_Seeding: TMainMenu
    Left = 208
    Top = 3
    object Refresh1: TMenuItem
      Bitmap.Data = {
        36030000424D3603000000000000360000002800000010000000100000000100
        18000000000000030000120B0000120B00000000000000000000FF00FFFF00FF
        FF00FFFF00FFFF00FFA23F08FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00
        FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFA13E08B04D06A13E08FF
        00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF
        FF00FFA13E08BD5804BF5B03A13E08FF00FFFF00FFFF00FFFF00FFFF00FFFF00
        FFFF00FFFF00FFFF00FFFF00FFFF00FFA13D08B34F05C86401A13E08A13E08FF
        00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFA13D08
        AE4C09C96502A44108A13E08FF00FFA13E08A13E08A03D08A03D08A03D08A03D
        08A03D08A13E08A13E08FF00FFA44209D97E1CAA4809A13E08FF00FFA13E08B0
        4D06C66102C66002C56002C56003C45F03C45F03C56003A03E08A03C06CE7626
        C96F1FA03C07FF00FFFF00FFFF00FFA13E08A13E08A13E08A13E08AD4B06CC68
        00C96402D06B00A03D089F3B06ED9D40A13E08FF00FFFF00FFFF00FFFF00FFFF
        00FFFF00FFFF00FF9F3C08BD5904C66202AA4606D06B00A03D089E3A05FFBF55
        9D3905FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFA03D08B24E05CD6900A23F
        09A13E08D06B00A03D089E3B05FAAE4BB151149F3C06FF00FFFF00FFFF00FFFF
        00FFA03D08BC5703CE6900A94607A03D08A13E08D16C00A03D089F3C06D27C2D
        E6953C9D3905A13E08FF00FF9F3C06A34109CC6D13CB6806A64207A13E08FF00
        FFA13E08D16C00A03D08FF00FFA5430AEE9E41E7953CB85A18B95A19C56A21ED
        9A37CA6E1BA44109A13E08FF00FFFF00FFA13E08D16C00A03D08FF00FFA03D07
        9F3B06BF631EE08E38EE9E41DA8532B15213A03C07A13E08FF00FFFF00FFFF00
        FFA13E08C96401A13E08FF00FFFF00FFFF00FFA03C069F3C069F3B069F3C06A0
        3D07FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFA13E08FF00FFFF00FFFF00FF
        FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00
        FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF
        00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF}
      Caption = 'Refresh'
      OnClick = Refresh1Click
    end
  end
  object SaveTextFileDialog1: TSaveTextFileDialog
    DefaultExt = 'txt'
    Left = 48
    Top = 372
  end
end
