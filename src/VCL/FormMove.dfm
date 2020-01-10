object frmMove: TfrmMove
  Left = 0
  Top = 0
  Width = 895
  Height = 637
  AutoScroll = True
  Caption = 'Movement Processes '
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -13
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 16
  object BtnOK: TButton
    Left = 680
    Top = 560
    Width = 89
    Height = 25
    Caption = 'OK'
    Default = True
    ModalResult = 1
    TabOrder = 5
    OnClick = BtnOKClick
  end
  object BtnCancel: TButton
    Left = 780
    Top = 560
    Width = 89
    Height = 25
    Cancel = True
    Caption = 'Cancel'
    ModalResult = 2
    TabOrder = 6
    OnClick = BtnCancelClick
  end
  object PanelRW: TPanel
    Left = 530
    Top = 328
    Width = 338
    Height = 218
    Alignment = taLeftJustify
    BevelInner = bvLowered
    Caption = 'Random Walks'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold, fsItalic]
    ParentFont = False
    TabOrder = 4
    VerticalAlignment = taAlignTop
    Visible = False
    object edtStepLength: TLabeledEdit
      Left = 109
      Top = 118
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 94
      EditLabel.Height = 16
      EditLabel.Caption = 'Step length (m) '
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
      Text = '0.0'
    end
    object edtRho: TLabeledEdit
      Left = 109
      Top = 148
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 67
      EditLabel.Height = 16
      EditLabel.Caption = 'Correlation '
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
      Text = '0.0'
    end
    object CBIndVarCRW: TCheckBox
      Left = 16
      Top = 33
      Width = 157
      Height = 17
      Caption = 'Individual variability'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold]
      ParentFont = False
      TabOrder = 0
      OnClick = CBIndVarCRWClick
    end
    object edtStepLMean: TLabeledEdit
      Left = 165
      Top = 118
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 36
      EditLabel.Height = 16
      EditLabel.Caption = 'Mean:'
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
      ParentFont = False
      TabOrder = 2
      Text = '0.0'
      Visible = False
    end
    object edtStepLSD: TLabeledEdit
      Left = 221
      Top = 118
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 28
      EditLabel.Height = 16
      EditLabel.Caption = 'S.d.:'
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
      ParentFont = False
      TabOrder = 3
      Text = '0.0'
      Visible = False
    end
    object edtStepMort: TLabeledEdit
      Left = 109
      Top = 178
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 84
      EditLabel.Height = 16
      EditLabel.Caption = 'Step mortality '
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
      TabOrder = 9
      Text = '0.0'
    end
    object edtRhoMean: TEdit
      Left = 165
      Top = 148
      Width = 50
      Height = 24
      Alignment = taRightJustify
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 6
      Text = '0.0'
      Visible = False
    end
    object edtRhoSD: TEdit
      Left = 221
      Top = 148
      Width = 50
      Height = 24
      Alignment = taRightJustify
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 7
      Text = '0.0'
      Visible = False
    end
    object edtStepLScale: TLabeledEdit
      Left = 277
      Top = 118
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.AlignWithMargins = True
      EditLabel.Width = 45
      EditLabel.Height = 32
      EditLabel.BiDiMode = bdLeftToRight
      EditLabel.Caption = 'Scaling factor:'
      EditLabel.Font.Charset = DEFAULT_CHARSET
      EditLabel.Font.Color = clWindowText
      EditLabel.Font.Height = -13
      EditLabel.Font.Name = 'Tahoma'
      EditLabel.Font.Style = []
      EditLabel.ParentBiDiMode = False
      EditLabel.ParentFont = False
      EditLabel.WordWrap = True
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 4
      Text = '0.0'
      Visible = False
    end
    object edtRhoScale: TEdit
      Left = 277
      Top = 148
      Width = 50
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
      Text = '0.0'
      Visible = False
    end
  end
  object PanelHab: TPanel
    Left = 530
    Top = 16
    Width = 337
    Height = 306
    Alignment = taLeftJustify
    BevelInner = bvLowered
    Caption = 'Habitat Costs / Mortality'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold, fsItalic]
    ParentFont = False
    TabOrder = 3
    VerticalAlignment = taAlignTop
    object SGhabCosts: TStringGrid
      Left = 8
      Top = 55
      Width = 198
      Height = 234
      BiDiMode = bdLeftToRight
      ColCount = 2
      DefaultColWidth = 60
      DoubleBuffered = False
      RowCount = 12
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -11
      Font.Name = 'Tahoma'
      Font.Style = []
      Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
      ParentBiDiMode = False
      ParentDoubleBuffered = False
      ParentFont = False
      TabOrder = 2
      ColWidths = (
        60
        60)
      RowHeights = (
        24
        24
        24
        24
        24
        24
        24
        24
        24
        24
        24
        24)
    end
    object CBCosts: TCheckBox
      Left = 8
      Top = 32
      Width = 169
      Height = 17
      Caption = 'Import cost map'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 0
      OnClick = CBCostsClick
    end
    object CBVisualCosts: TCheckBox
      Left = 152
      Top = 33
      Width = 195
      Height = 16
      Caption = 'Visualise costs landscape'
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 1
      Visible = False
      WordWrap = True
    end
  end
  object RGMoveModel: TRadioGroup
    Left = 8
    Top = 8
    Width = 249
    Height = 81
    Caption = 'Movement Model'
    DoubleBuffered = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ItemIndex = 0
    Items.Strings = (
      'SMS - Stochastic Movement Simulator'
      'CRW (continuous space)')
    ParentDoubleBuffered = False
    ParentFont = False
    TabOrder = 0
    WordWrap = True
    OnClick = RGMoveModelClick
  end
  object RGstepM: TRadioGroup
    Left = 8
    Top = 104
    Width = 177
    Height = 73
    Caption = 'Step Mortality'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ItemIndex = 0
    Items.Strings = (
      'Constant'
      'Habitat dependent')
    ParentFont = False
    TabOrder = 1
    OnClick = RGstepMClick
  end
  object PanelSMS: TPanel
    Left = 8
    Top = 183
    Width = 505
    Height = 363
    Alignment = taLeftJustify
    BevelInner = bvLowered
    Caption = 'SMS'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold, fsItalic]
    ParentFont = False
    TabOrder = 2
    VerticalAlignment = taAlignTop
    object PRLabel: TLabel
      Left = 205
      Top = 52
      Width = 68
      Height = 16
      Caption = '>= 1 (cells)'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
    end
    object PRmethodLabel: TLabel
      Left = 205
      Top = 82
      Width = 59
      Height = 16
      Caption = '1 / 2 / 3 *'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
    end
    object DPLabel: TLabel
      Left = 205
      Top = 113
      Width = 44
      Height = 16
      Caption = '>= 1.0 '
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
    end
    object GoalBiasLabel: TLabel
      Left = 205
      Top = 172
      Width = 44
      Height = 16
      Caption = '>= 1.0 '
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
    end
    object MemorySizeLabel: TLabel
      Left = 205
      Top = 141
      Width = 34
      Height = 16
      Caption = '1 - 14'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
    end
    object GoalTypeLabel: TLabel
      Left = 205
      Top = 201
      Width = 47
      Height = 16
      Caption = '0 / 1 / 2'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
    end
    object AlphaDBlabel: TLabel
      Left = 205
      Top = 231
      Width = 35
      Height = 16
      Caption = '> 0.0 '
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
    end
    object BetaDBlabel: TLabel
      Left = 205
      Top = 261
      Width = 68
      Height = 16
      Caption = '>= 1 (cells)'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
    end
    object edtPRMethod: TLabeledEdit
      Left = 149
      Top = 78
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 66
      EditLabel.Height = 16
      EditLabel.Caption = 'PR method '
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
      Text = '1'
    end
    object edtPR: TLabeledEdit
      Left = 149
      Top = 48
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 140
      EditLabel.Height = 16
      EditLabel.Caption = 'Perceptual range (cells) '
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
      TabOrder = 11
      Text = '1'
    end
    object edtDP: TLabeledEdit
      Left = 149
      Top = 108
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 133
      EditLabel.Height = 16
      EditLabel.Caption = 'Directional persistence '
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
      Text = '5.0'
    end
    object edtSM: TLabeledEdit
      Left = 149
      Top = 18
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 84
      EditLabel.Height = 16
      EditLabel.Caption = 'Step mortality '
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
      TabOrder = 12
      Text = '0.0'
    end
    object Memo2: TMemo
      Left = 80
      Top = 300
      Width = 200
      Height = 63
      BevelKind = bkFlat
      BorderStyle = bsNone
      Color = clInfoBk
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsItalic]
      Lines.Strings = (
        '* 1 = arithmetic mean, '
        '   2 = harmonic mean, '
        '   3 = weighted arithmetic mean')
      ParentFont = False
      ReadOnly = True
      TabOrder = 7
    end
    object edtMemorySize: TLabeledEdit
      Left = 149
      Top = 138
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 76
      EditLabel.Height = 16
      EditLabel.Caption = 'Memory size '
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
      TabOrder = 4
      Text = '1'
    end
    object edtGoalBias: TLabeledEdit
      Left = 149
      Top = 168
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 56
      EditLabel.Height = 16
      EditLabel.Caption = 'Goal bias '
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
      Text = '1.0'
    end
    object edtGoalType: TLabeledEdit
      Left = 149
      Top = 198
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 57
      EditLabel.Height = 16
      EditLabel.Caption = 'Goal type '
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
      Text = '0'
      OnChange = edtGoalTypeChange
    end
    object edtAlphaDB: TLabeledEdit
      Left = 149
      Top = 228
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 55
      EditLabel.Height = 16
      EditLabel.Caption = 'Alpha DB '
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
      Text = '1.0'
    end
    object edtBetaDB: TLabeledEdit
      Left = 149
      Top = 258
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 48
      EditLabel.Height = 16
      EditLabel.Caption = 'Beta DB '
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
      Text = '100000'
    end
    object CBIndVarSMS: TCheckBox
      Left = 300
      Top = 22
      Width = 161
      Height = 17
      Caption = 'Individual variability'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold]
      ParentFont = False
      TabOrder = 9
      OnClick = CBIndVarSMSClick
    end
    object PanelDBIndVar: TPanel
      Left = 290
      Top = 224
      Width = 200
      Height = 70
      BevelOuter = bvNone
      Enabled = False
      TabOrder = 10
      Visible = False
      object edtAlphaDBMean: TEdit
        Left = 10
        Top = 4
        Width = 50
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 0
        Text = '0.0'
      end
      object edtAlphaDBSD: TEdit
        Left = 70
        Top = 4
        Width = 50
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 1
        Text = '0.0'
      end
      object edtAlphaDBScale: TEdit
        Left = 130
        Top = 4
        Width = 50
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 2
        Text = '0.0'
      end
      object edtBetaDBMean: TEdit
        Left = 10
        Top = 34
        Width = 50
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 3
        Text = '0'
      end
      object edtBetaDBSD: TEdit
        Left = 70
        Top = 34
        Width = 50
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 4
        Text = '0'
      end
      object edtBetaDBScale: TEdit
        Left = 130
        Top = 34
        Width = 50
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 5
        Text = '0'
      end
    end
    object PanelDPGBIndVar: TPanel
      Left = 290
      Top = 74
      Width = 200
      Height = 122
      BevelOuter = bvNone
      Enabled = False
      TabOrder = 8
      Visible = False
      object edtGBMean: TEdit
        Left = 10
        Top = 95
        Width = 50
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 3
        Text = '0.0'
      end
      object edtGBSD: TEdit
        Left = 70
        Top = 95
        Width = 50
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 4
        Text = '0.0'
      end
      object edtGBScale: TEdit
        Left = 130
        Top = 95
        Width = 50
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 5
        Text = '0.0'
      end
      object edtDPMean: TLabeledEdit
        Left = 10
        Top = 35
        Width = 50
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 36
        EditLabel.Height = 16
        EditLabel.Caption = 'Mean:'
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
        ParentFont = False
        TabOrder = 0
        Text = '0.0'
      end
      object edtDPSD: TLabeledEdit
        Left = 70
        Top = 35
        Width = 50
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 28
        EditLabel.Height = 16
        EditLabel.Caption = 'S.d.:'
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
        ParentFont = False
        TabOrder = 1
        Text = '0.0'
      end
      object edtDPScale: TLabeledEdit
        Left = 130
        Top = 35
        Width = 50
        Height = 24
        Alignment = taRightJustify
        EditLabel.AlignWithMargins = True
        EditLabel.Width = 45
        EditLabel.Height = 32
        EditLabel.BiDiMode = bdLeftToRight
        EditLabel.Caption = 'Scaling factor:'
        EditLabel.Font.Charset = DEFAULT_CHARSET
        EditLabel.Font.Color = clWindowText
        EditLabel.Font.Height = -13
        EditLabel.Font.Name = 'Tahoma'
        EditLabel.Font.Style = []
        EditLabel.ParentBiDiMode = False
        EditLabel.ParentFont = False
        EditLabel.WordWrap = True
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 2
        Text = '0.0'
      end
    end
  end
  object CBstraightPath: TCheckBox
    Left = 8
    Top = 560
    Width = 281
    Height = 23
    Caption = 'Straighten path after decision not to settle'
    Checked = True
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentFont = False
    State = cbChecked
    TabOrder = 7
  end
  object OpenDialog1: TOpenDialog
    Left = 240
    Top = 112
  end
end
