define([
    'app/app',
    'text!app/templates/user/list.hbs',
    'app/view/user/thumb'
], function (App, template) {

    var $window = $(window);

    App.View.User = App.View.User || {};

    App.View.User.List = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thread-members thread-box',

        invite_button: false,
        template: App.Handlebars.compile(template),

        $list: null,
        resizeTimeout: null,
        nbThumbsPerRows: 1,
        expandedNbRows: 1,

        events: {
            'click .thread-members-toggle a': 'onToggleClick'
        },

        initialize: function (options) {
            App._.defaults(options, {invite_button: false});

            this.invite_button = options.invite_button;

            this.listenTo(this.model, "reset update", this.render);

            App._.bindAll(this, 'onWindowResize', 'delayOnWindowResize');
            $window.on('resize', this.delayOnWindowResize);
        },

        remove: function() {
            $window.off('resize', this.delayOnWindowResize);

            App.Backbone.View.prototype.remove.apply(this, arguments);
        },

        _calculateNbThumbsPerRow: function() {
            // thumb = 84x84
            return Math.floor(this.$el.width() / 84);
        },

        _calculateNbRows: function() {
            return Math.ceil(this.model.size() / this._calculateNbThumbsPerRow());
        },

        delayOnWindowResize: function() {
            if (this.resizeTimeout) {
                clearTimeout(this.resizeTimeout);
            }
            this.resizeTimeout = setTimeout(this.onWindowResize, 250);
        },

        onWindowResize: function() {
            this.resizeTimeout = null;

            this.nbThumbsPerRows = this._calculateNbThumbsPerRow();
            this.assignThumbsDirection();

            this.expandedNbRows = this._calculateNbRows();
            this.assignListHeight();
        },

        assignThumbsDirection: function() {
            if (0 == this.nbThumbsPerRows) {
                return;
            }
            var $thumbs = this.$('.user').removeClass('alt');
            $thumbs.filter(':nth-child(' + this.nbThumbsPerRows + 'n)').addClass('alt');
            $thumbs.filter('.user:nth-child(' + this.nbThumbsPerRows + 'n+' + (this.nbThumbsPerRows - 1) + ')').addClass('alt');
        },

        assignListHeight: function() {
            var $icon = this.$('.thread-members-toggle a span');

            if (1 < this.expandedNbRows) {
                $icon.show();
                if (this.$list.hasClass('expanded')) {
                    this.$list.css({height: (this.expandedNbRows * 84)});
                    $icon.removeClass('eai-down').addClass('eai-up');
                    return;
                }
            } else {
                $icon.hide();
            }

            this.$list.removeClass('expanded').css({height: 80});
            $icon.removeClass('eai-up').addClass('eai-down');
        },

        onToggleClick: function(e) {
            e.preventDefault();

            if (1 < this.expandedNbRows && !this.$list.hasClass('expanded')) {
                this.$list.addClass('expanded')
            } else {
                this.$list.removeClass('expanded')
            }
            this.assignListHeight();
        },

        render: function () {
            App.log('UserListView::render');

            this.$el.html(this.template({
                invite_button: this.invite_button
            }));

            var that = this;
            this.$list = this.$('.thread-members-list');

            this.model.each(function(user) {
                that.pushSubView(
                    App.View.User.Thumb.create({
                        model: user
                    }, {
                        $target: that.$list,
                        append: true
                    })
                );
            });

            this.delayOnWindowResize();

            return this;
        }
    });

    App.View.User.List.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            options.model = new App.Model.Users();
        }

        var list = new App.View.User.List(options);
        if (createOptions) {
            App.View.create(list, createOptions);
        }

        return list;
    };

    return App.View.User.List;
});
