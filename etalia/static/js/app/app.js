define([
    'require',
    'jquery',
    'underscore',
    'backbone',
    'handlebars',
    'moment',
    'app/ui/layout',
    'bootstrap',
    'backbone/relational',
    'backbone/paginator',
    'backbone/forms'
], function (require, $, _, Backbone, Handlebars, Moment, Layout) {

    /**
     * Backbone
     */
    var _sync = Backbone.sync;
    Backbone.sync = function (method, model, options) {
        // Add trailing slash to backbone model views
        var _url = options.url ? options.url : _.isFunction(model.url) ? model.url() : model.url;
        _url += _url.charAt(_url.length - 1) == '/' ? '' : '/';
        options = _.extend(options, {
            url: _url
        });
        return _sync(method, model, options);
    };

    var BBColFetch = Backbone.Collection.prototype.fetch;
    Backbone.Collection.prototype.fetch = function (options) {
        options = options || {};
        if (!options.data && this.queryParams) {
            options.data = this.queryParams;
        }
        return BBColFetch.call(this, options);
    };
    Backbone.Collection.prototype.parse = function (response, options) {
        if (_.has(response, 'results')) {
            return response.results;
        }
        return response;
    };

    var BBViewRemove = Backbone.View.prototype.remove;
    Backbone.View.prototype.pushSubView = function (view) {
        if (!this.subViews) {
            this.subViews = [];
        }
        this.subViews.push(view);
    };
    Backbone.View.prototype.clearSubViews = function() {
        if (this.subViews) {
            _.each(this.subViews, function (view) {
                view.remove();
            });
            this.subViews = null;
        }
    };
    Backbone.View.prototype.remove = function () {
        this.clearSubViews();
        BBViewRemove.apply(this, arguments);
    };

    /*Backbone.View.prototype.__remove = Backbone.View.prototype.remove;
    Backbone.View = Backbone.View.extend({
        pushSubView: function (view) {
            if (!this.subViews) {
                this.subViews = [];
            }
            this.subViews.push(view);
        },
        clearSubViews: function() {
            if (this.subViews) {
                _.each(this.subViews, function (view) {
                    view.remove();
                });
                this.subViews = null;
            }
        },
        remove: function () {
            this.clearSubViews();
            Backbone.View.prototype.__remove.apply(this, arguments);
        }
    });*/

    /**
     * Backbone Relational
     */
    function isApiUrl(url) {
        return _.isString(url) && /https?:\/\/[0-9a-z-:\.]+\/[a-z0-9-\/]+\/[0-9]+\//.test(url);
    }

    function extractIdFromApiUrl(url) {
        // TODO assert that the ID is properly extracted.
        return parseInt(/\d+/.exec(/\/\d+\//.exec(url)[0])[0]);
    }

    var relationalModel_updateRelations =
        Backbone.RelationalModel.prototype.__updateRelations =
            Backbone.RelationalModel.prototype.updateRelations;
    Backbone.RelationalModel.prototype.updateRelations = function (changedAttrs, options) {
        if (this._isInitialized && !this.isLocked()) {
            _.each(this._relations, function (rel) {
                var value;
                if (rel.keySource in changedAttrs) {
                    value = this.attributes[rel.keySource];
                    if (isApiUrl(value)) {
                        this.attributes[rel.keySource] = extractIdFromApiUrl(value);
                        changedAttrs[rel.keySource] = extractIdFromApiUrl(value);
                    }
                }
                if (rel.key in changedAttrs) {
                    value = this.attributes[rel.key];
                    if (isApiUrl(value)) {
                        this.attributes[rel.key] = extractIdFromApiUrl(value);
                        changedAttrs[rel.key] = extractIdFromApiUrl(value);
                    }
                }
            }, this);
        }

        relationalModel_updateRelations.call(this, changedAttrs, options);
    };

    var hasOneSetKeyContents = Backbone.HasOne.prototype.__setKeyContents = Backbone.HasOne.prototype.setKeyContents;
    Backbone.HasOne.prototype.setKeyContents = function (keyContents) {
        // Extract ID from Api Url
        if (isApiUrl(keyContents)) {
            keyContents = extractIdFromApiUrl(keyContents);
        }
        hasOneSetKeyContents.call(this, keyContents);
    };

    // @TODO test
    var hasManySetKeyContents = Backbone.HasMany.prototype.__setKeyContents = Backbone.HasMany.prototype.setKeyContents;
    Backbone.HasMany.prototype.setKeyContents = function (keyContents) {

        // Asserts that keyContents is an array of api url string
        if (_.isArray(keyContents)) {
            var assert = true;
            _.each(keyContents, function (value) {
                if (!isApiUrl(value)) {
                    assert = false;
                }
            });
        }

        // Convert url[] to id[]
        if (assert) {
            keyContents = _.map(keyContents, function (value) {
                return extractIdFromApiUrl(value);
            });
        }

        hasManySetKeyContents.call(this, keyContents);
    };
    Backbone.HasMany.prototype.getCount = function () {
        var relatedCount = this.related.models.length,
            idsCount = this.keyIds.length;

        return relatedCount > 0 ? relatedCount : idsCount;
    };


    /**
     * Backbone Paginator
     */
    Backbone.PageableCollection.prototype.parseRecords = function (response) {
        if (_.has(response, 'results')) {
            return response.results;
        }
        return response;
    };
    Backbone.PageableCollection.prototype.parseLinks = function (response) {
        var links = {};
        if (_.has(response, 'next')) {
            links.next = response.next;
        }
        if (_.has(response, 'prev')) {
            links.prev = response.prev;
        }
        return links;
    };


    /**
     * Backbone Forms
     */
        // @see https://github.com/gerardobort/bbf-bbr-integration/blob/master/test.js
    Backbone.Form.editors.Base.prototype._initialize = Backbone.Form.editors.Base.prototype.initialize;
    Backbone.Form.editors.Base.prototype.initialize = function (options) {
        Backbone.Form.editors.Base.prototype._initialize.call(this, options);
        // this patch adds compatibility between backbone-forms and backbone-relational libraries
        if (options.model instanceof Backbone.RelationalModel && options.model.get(options.key) instanceof Backbone.Collection) {
            this.value = options.model.get(options.key).toJSON();
        }
    };

    Backbone.Form.template = _.template('\
    <form class="form-horizontal" role="form">\
      <div data-fieldsets></div>\
      <% if (submitButton) { %>\
        <button type="submit" class="btn btn-primary"><%= submitButton %></button>\
      <% } %>\
    </form>\
    ');
    Backbone.Form.Field.template = _.template('\
    <div class="form-group field-<%= key %>">\
      <% if (title) { %>\
      <label class="col-sm-2 control-label hidden-xs" for="<%= editorId %>"><%= title %></label>\
      <% } %>\
      <div class="col-sm-<% if (title) { %>10<% } else { %>12<% } %>">\
        <span data-editor></span>\
        <p class="help-block" data-error></p>\
        <p class="help-block"><%= help %></p>\
      </div>\
    </div>\
    ');

    Backbone.Form.editors.Select.prototype.__setValue = Backbone.Form.editors.Select.prototype.setValue;
    Backbone.Form.editors.Select = Backbone.Form.editors.Select.extend({
        setValue: function (value) {
            if (typeof value === 'object') {
                try {
                    value = value.get('id');
                } catch (e) {
                }
            }
            this.__setValue(value);
        }
    });

    Backbone.Form.editors.Radio.prototype.__setValue = Backbone.Form.editors.Radio.prototype.setValue;
    Backbone.Form.editors.Radio = Backbone.Form.editors.Radio.extend({
        tagName: 'div',
        className: 'btn-group',
        attributes: {
            'data-toggle': 'buttons'
        },
        setValue: function (value) {
            this.__setValue(value);
            this.$('input[type=radio]:checked').closest('label').addClass('active');
        }
    }, {
        template: _.template('\
    <% _.each(items, function(item) { %>\
      <label for="<%= item.id %>" class="btn btn-default">\
        <input type="radio" name="<%= item.name %>" value="<%= item.value %>" id="<%= item.id %>" autocomplete="off"> <%= item.label %>\
      </label>\
    <% }); %>\
    ', null, Backbone.Form.templateSettings)
    });

    Backbone.Form.editors.Tinymce = Backbone.Form.editors.TextArea.extend({

        /*tagName: 'textarea',

         events: {
         'change': function() {
         // The 'change' event should be triggered whenever something happens
         // that affects the result of `this.getValue()`.
         this.trigger('change', this);
         },
         'focus': function() {
         // The 'focus' event should be triggered whenever an input within
         // this editor becomes the `document.activeElement`.
         this.trigger('focus', this);
         // This call automatically sets `this.hasFocus` to `true`.
         },
         'blur': function() {
         // The 'blur' event should be triggered whenever an input within
         // this editor stops being the `document.activeElement`.
         this.trigger('blur', this);
         // This call automatically sets `this.hasFocus` to `false`.
         }
         },*/

        defaultSettings: {
            resize: false,
            schema: "html5-strict",
            fix_list_elements: true,
            keep_styles: false,
            valid_elements: "h1,h2,h3,h4,h5,h6,p,br,hr,ul,ol,li,blockquote,pre,strong/b,em/i,a[!href|target=_blank]",
            statusbar: false,
            menubar: false,
            height: 160,
            plugins: "lists,paste,link,searchreplace,hr",
            toolbar: "undo redo | formatselect | bold italic | bullist numlist | link hr blockquote",
            skin_url: '/static/css/lib/tinymce',
            body_class: 'element-content',
            content_css: [
                'https://fonts.googleapis.com/css?family=Lato:400,700,700italic,400italic',
                '/static/css/lib/bootstrap.css',
                '/static/css/app/content.css'
            ]
        },

        tinymceEditorId: null,

        initialize: function (options) {
            Backbone.Form.editors.Base.prototype.initialize.call(this, options);
        },

        render: function () {
            this.setValue(this.value);

            var that = this,
                settings = {};

            this.tinymceEditorId = that.$el.attr('id');

            _.defaults(
                settings,
                this.schema.settings,
                this.defaultSettings,
                {selector: '#' + that.tinymceEditorId}
            );
            settings.setup = function (editor) {
                var onChange = function () {
                    that.setValue(editor.getContent());

                    /*und.convert(editor.getContent(), function(err, markdown) {
                     $('#markdown').val(markdown);
                     }, { keepHtml: false });*/
                };

                editor.on('keyup', onChange);
                editor.on('change', onChange);
            };

            require(['tinymce'], function () {
                tinymce.EditorManager.execCommand('mceRemoveEditor', true, that.tinymceEditorId);

                tinymce.init(settings);
            });

            return this;
        },

        remove: function () {
            var that = this;
            require(['tinymce'], function () {
                tinymce.EditorManager.execCommand('mceRemoveEditor', true, that.tinymceEditorId);
            });

            if (this.tinymceEditor) {
                this.tinymceEditor.destroy();
                this.tinymceEditor = null;
            }

            Backbone.View.prototype.remove.call(this);
        }

        /*getValue: function() {
         return this.$el.val();
         },

         setValue: function(value) {
         this.$el.val(value);
         },

         focus: function() {
         if (this.hasFocus) return;

         // This method call should result in an input within this editor
         // becoming the `document.activeElement`.
         // This, in turn, should result in this editor's `focus` event
         // being triggered, setting `this.hasFocus` to `true`.
         // See above for more detail.
         this.$el.focus();
         },

         blur: function() {
         if (!this.hasFocus) return;

         this.$el.blur();
         }*/
    });


    /**
     * Handlebars helpers
     */

    Handlebars.registerHelper('ifCond', function (a, operator, b, options) {
        switch (operator) {
            case '==':
                return (a == b) ? options.fn(this) : options.inverse(this);
            case '===':
                return (a === b) ? options.fn(this) : options.inverse(this);
            case '<':
                return (a < b) ? options.fn(this) : options.inverse(this);
            case '<=':
                return (a <= b) ? options.fn(this) : options.inverse(this);
            case '>':
                return (a > b) ? options.fn(this) : options.inverse(this);
            case '>=':
                return (a >= b) ? options.fn(this) : options.inverse(this);
            case 'and':
                return (a && b) ? options.fn(this) : options.inverse(this);
            case 'andNot':
                return (a && !b) ? options.fn(this) : options.inverse(this);
            case 'or':
                return (a || b) ? options.fn(this) : options.inverse(this);
            case 'orNot':
                return (a || !b) ? options.fn(this) : options.inverse(this);
            default:
                return options.inverse(this);
        }
    });
    Handlebars.registerHelper('formatDate', function (date, format) {
        if (!_.isString(format)) {
            format = 'MMM D, YYYY';
        }
        return Moment(date).format(format);
    });
    Handlebars.registerHelper('loadingSprite', function () {
        return new Handlebars.SafeString(
            '<div class="sk-cube-grid">'+
                '<div class="sk-cube sk-cube1"></div>'+
                '<div class="sk-cube sk-cube2"></div>'+
                '<div class="sk-cube sk-cube3"></div>'+
                '<div class="sk-cube sk-cube4"></div>'+
                '<div class="sk-cube sk-cube5"></div>'+
                '<div class="sk-cube sk-cube6"></div>'+
                '<div class="sk-cube sk-cube7"></div>'+
                '<div class="sk-cube sk-cube8"></div>'+
                '<div class="sk-cube sk-cube9"></div>'+
            '</div>'
        );
    });


    /**
     * Polyfills
     */
    if (!String.prototype.trim) {
        String.prototype.trim = function (str) {
            return str.replace(/^\s+|\s+$/gm,'');
        }
    }

    var Validator = {
        isRelation: function(relation) {
            if (_.isObject(relation)) {
                if (relation instanceof Backbone.Model) {
                    relation = relation.attributes;
                }
                if (relation.hasOwnProperty('link') && 0 < String(relation.link).length) {
                    return true;
                }
            } else if (0 < parseInt(relation)) {
                return true;
            }
            return false;
        }
    };

    /*$.fn.isFixed = function() {
        if (1 !== this.size()) {
            throw 'Collection has more or less than one element.';
        }
        if (this.css('position') === 'fixed') return true;
        var fixed = false;
        this.parents().each(function(){
            if ($(this).css('position') === 'fixed') {
                fixed = true;
                return false;
            }
        });
        return fixed;
    };*/

    /**
     * App
     */
    window.App = App = _.extend({
        $: $,
        _: _,
        Backbone: Backbone,
        Handlebars: Handlebars,

        Validator: Validator,
        Layout: Layout,

        Const: {},
        Model: {},
        Collection: {},
        View: {
            createDefaults: {
                $target: null,
                append: false
            },
            create: function (view, options) {
                options = _.defaults(options, App.View.createDefaults);
                if (options.$target && options.$target.length) {
                    if (true === options.append) {
                        options.$target.append(view.render().$el);
                    } else {
                        options.$target.replaceWith(view.render().$el);
                    }
                }
            }
        },

        config: {
            debug: false,
            api_root: '/api/v1'
        },

        defaultHeadTitle: 'Etalia',

        log: function () {
            if (this.config.debug) {
                if (arguments.length > 1) {
                    console.log(arguments[0], Array.prototype.splice.call(arguments, 1));
                } else {
                    console.log(arguments[0]);
                }
            }
        },

        popup: function(url, name, width, height) {
            name = name || 'etalia-popup';
            width = width || 520;
            height = height || 460;

            var left = Math.round((window.innerWidth/ 2) - (width / 2)),
                top = Math.round((window.innerHeight / 2) - (height / 2)),
                params = "menubar=no,toolbar=no,resizable=yes,scrollbars=yes," +
                    "width=" + width + ",height=" + height + "," +
                    "top=" + top + ",left=" + left;

            return window.open(url, name, params);
        },

        getProperty: function(object, path) {
            return object.hasOwnProperty(path) ? object[path] : object.get(path);
        }
    }, Backbone.Events);

    App.loadPopovers = function() {
        // Close all on click out
        $(document).on('click', function(e) {
            if ($(e.target).closest('.ui-user-popover').size() == 0 && $(e.target).closest('.popover').size() == 0) {
                $('.ui-user-popover').popover('hide');
            }
        });

        require(['app/view/ui/user-popover', 'app/view/ui/modal'], function() {

            var userPopovers = new App.Model.UserPopovers();
            userPopovers.on('change', _loadPopovers);
            _loadPopovers();

            function _loadPopovers() {
                userPopovers.fetch().then(function() {
                    var modalUserPopover = userPopovers.find(function (userPopover) {
                        return userPopover.get('popover').get('type') === App.Model.Popover.TYPE_MODAL;
                    });
                    if (modalUserPopover) {
                        var data = modalUserPopover.get('popover');
                        var modal = new App.View.Ui.Modal({
                            title: data.get('title'),
                            content: data.get('body'),
                            footer: {
                                buttons: [{
                                    label: 'Got it !',
                                    attr: {
                                        'type': 'button',
                                        'class': 'btn btn-primary popover-got-it'
                                    }
                                }]
                            }
                        });
                        modal
                            .render()
                            .on('button_click', function (e) {
                                var $button = $(e.target).closest('button');
                                if ($button.hasClass('popover-got-it')) {
                                    modalUserPopover.markDone().then(function() {
                                        modal.close();
                                    });
                                }
                            });
                    } else {
                        userPopovers.each(function (userPopover) {
                            if (userPopover.get('popover').get('type') === App.Model.Popover.TYPE_ANCHORED) {
                                var $target = $(userPopover.get('popover').get('anchor')).first();
                                if (1 == $target.size() && $target.is(':visible')) {
                                    var userPopoverView = App.View.Ui.UserPopover.create({
                                        model: userPopover,
                                        $target: $target
                                    });
                                }
                            } else {
                                throw 'Unexpected popover type';
                            }
                        });
                    }
                });
            }
        });
    };

    App.loadTracking = function() {
        require(['app/util/tracking'], function(Tracking) {
            Tracking.init();
        });
    };

    App.defaultHeadTitle = $('title').text();

    App.setHeadTitle = function(title) {
        $('title').text(!!title ? title : App.defaultHeadTitle);
    };

    App.init = function() {
        App.loadPopovers();
        App.loadTracking();
    };

    return App;
});
