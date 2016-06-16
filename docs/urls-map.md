URL maps
========

/feed/stream/	etalia.feeds.views.StreamView	feeds:home	
/feed/stream/	etalia.feeds.views.StreamView	feeds:stream	
/feed/stream/<name>/reset	etalia.feeds.views.reset_stream_view	feeds:reset-stream	login_required
/feed/stream/<name>/update	etalia.feeds.views.update_stream_view	feeds:update-stream	login_required
/feed/stream/xml	etalia.feeds.views.StreamViewXML	feeds:stream_xml	
/feed/trend/	etalia.feeds.views.TrendView	feeds:trend	
/feed/trend/<name>/reset	etalia.feeds.views.reset_trend_view	feeds:reset-trend	login_required
/feed/trend/<name>/update	etalia.feeds.views.update_trend_view	feeds:update-trend	login_required
/feed/trend/xml	etalia.feeds.views.TrendViewXML	feeds:trend_xml	
/help/	etalia.core.views.help	core:help	
/library/	etalia.library.views.library	library:library	
/library/journal/<pk>/	etalia.library.views.JournalViewPk	library:journal	
/library/journal/<slug>-<pk>/	etalia.library.views.JournalView	library:journal-slug	
/library/journals/	etalia.library.views.JournalsListView	library:journals	
/library/paper/<pk>/	etalia.library.views.PaperViewPk	library:paper	
/library/paper/<pk>/neighbors	etalia.library.views.PaperNeighborsView	library:paper-neighbors	
/library/paper/<slug>-<pk>/	etalia.library.views.PaperView	library:paper-slug	
/messages/mark_read/<message_id>/	messages_extends.views.message_mark_read	message_extends:message_mark_read	
/messages/mark_read/all/	messages_extends.views.message_mark_all_read	message_extends:message_mark_all_read	
/support/	etalia.core.views.support	core:support	
/terms/	etalia.core.views.terms	core:terms	
/test-failing-task	etalia.core.views.test_failing_task		
/threads/<pk>/	etalia.threads.views.ThreadView	threads:thread	
/threads/<pk>/join	etalia.threads.views.join_thread	threads:join_thread	
/threads/<pk>/posts/new	etalia.threads.views.ThreadPostCreateView	threads:new_post	
/threads/<pk>/update	etalia.threads.views.ThreadUpdateView	threads:update_thread	
/threads/comments/<pk>/delete	etalia.threads.views.ThreadPostCommentDeleteView	threads:delete_comment	
/threads/comments/<pk>/update	etalia.threads.views.ThreadPostCommentUpdateView	threads:update_comment	
/threads/new	etalia.threads.views.ThreadCreate	threads:new_thread	
/threads/posts/<pk>/comments/new	etalia.threads.views.ThreadPostCommentCreateView	threads:new_comment	
/threads/posts/<pk>/delete	etalia.threads.views.ThreadPostDeleteView	threads:delete_post	
/threads/posts/<pk>/update	etalia.threads.views.ThreadPostUpdateView	threads:update_post	
/user/avatar/add/	avatar.views.add	user:avatar_add	login_required
/user/avatar/change/	avatar.views.change	user:avatar_change	login_required
/user/avatar/delete/	avatar.views.delete	user:avatar_delete	login_required
/user/avatar/list/<username>/	avatar.views.avatar_gallery	user:avatar_gallery	login_required
/user/avatar/list/<username>/<id>/	avatar.views.avatar	user:avatar	login_required
/user/avatar/render_primary/<user>/<size>/	avatar.views.render_primary	user:avatar_render_primary	login_required
/user/complete/<backend>/	social.apps.django_app.views.complete	social:complete	
/user/delete	etalia.users.views.UserDeleteView	user:delete-user	
/user/disconnect/<backend>/	social.apps.django_app.views.disconnect	social:disconnect	
/user/disconnect/<backend>/<association_id>/	social.apps.django_app.views.disconnect	social:disconnect_individual	
/user/done/	etalia.users.views.done	user:done	login_required
/user/email-sent/	etalia.users.views.validation_sent	user:validation-sent	login_required
/user/info-affiliation/	etalia.users.views.UserAffiliationSignupView	user:require-affiliation	
/user/info/	etalia.users.views.UserBasicInfoSignupView	user:require-basic-info	
/user/library/	etalia.users.views.UserLibraryView	user:library	
/user/library/pins/	etalia.users.views.UserLibraryPinsView	user:library-pins	
/user/library/pins/xml	etalia.users.views.UserLibraryPinsViewXML	user:library-pins-xml	
/user/library/trash/	etalia.users.views.UserLibraryTrashView	user:library-trash	
/user/library/trash/empty	etalia.users.views.empty_trash_call	user:empty-trash	login_required
/user/library/trash/xml	etalia.users.views.UserLibraryTrashViewXML	user:library-trash-xml	
/user/library/xml	etalia.users.views.UserLibraryViewXML	user:library-xml	
/user/login/<backend>/	social.apps.django_app.views.auth	social:begin	
/user/logout/	etalia.users.views.logout	user:logout	login_required
/user/paper/add	etalia.users.views.AddCallView	user:add	
/user/paper/ban	etalia.users.views.BanCallView	user:ban	
/user/paper/pin	etalia.users.views.PinCallView	user:pin	
/user/paper/restore	etalia.users.views.RestoreCallView	user:restore	
/user/paper/trash	etalia.users.views.TrashCallView	user:trash	
/user/profile/	etalia.users.views.ProfileView	user:profile	
/user/send-invite	etalia.users.views.send_invite	user:send-invite	login_required
/user/settings/	etalia.users.views.SettingsView	user:settings	
/user/signin/	etalia.users.views.UserLoginView	user:signin	
/user/threads/	etalia.users.views.ThreadsView	user:threads	
/user/update-affiliation	etalia.users.views.UserAffiliationUpdateView	user:update-affiliation	
/user/update-email-digest-settings	etalia.users.views.UserEmailDigestSettingsUpdateView	user:update-email-digest-settings	
/user/update-library	etalia.users.views.update_library	user:update-library	login_required
/user/update-name	etalia.users.views.UpdateUserNameView	user:update-name	
/user/update-position	etalia.users.views.UpdateUserPositionView	user:update-position	
/user/update-stream-settings	etalia.users.views.UserStreamSettingsUpdateView	user:update-stream-settings	
/user/update-title	etalia.users.views.UpdateUserTitleView	user:update-title	
/user/update-trend-settings	etalia.users.views.UserTrendSettingsUpdateView	user:update-trend-settings	
/user/user-init-check	etalia.users.views.user_init_check	user:user-init-check	login_required
/user/user-update-settings-check	etalia.users.views.user_update_settings_check	user:user-update-settings-check	login_required
/user/user-update-stream-check	etalia.users.views.user_update_stream_check	user:user-update-stream-check	login_required
/user/user-update-trend-check	etalia.users.views.user_update_trend_check	user:user-update-trend-check	login_required

